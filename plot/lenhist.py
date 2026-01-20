#!/usr/bin/env -S uv run
# /// script
# dependencies = [
#   "pandas",
#   "altair",
#   "vl-convert-python",
# ]
# ///

import argparse
import json
import os
import sys

import altair as alt
import pandas as pd


def load_data(input_file, target=None):
    """Load length histogram data from JSON file.

    Args:
        input_file: Path to JSON file from 'grate len'
        target: Which histogram to use:
            - None or "total": Use total_stats.length_histogram (deduplicated)
            - "all": Use all per-target histograms
            - <name>: Use histogram for specific target

    Returns:
        DataFrame with columns: sample, length, count, [target]
    """
    with open(input_file, 'r') as f:
        data = json.load(f)

    rows = []

    for sample in data.get("samples", []):
        sample_name = sample.get("sample_name", "unknown")

        if target is None or target == "total":
            # Use total_stats histogram (deduplicated across targets)
            total_stats = sample.get("total_stats", {})
            length_histogram = total_stats.get("length_histogram", [])

            for length, count in length_histogram:
                rows.append({
                    "sample": sample_name,
                    "target": "Total",
                    "length": length,
                    "count": count
                })
        elif target == "all":
            # Use all per-target histograms
            for tgt in sample.get("targets", []):
                target_name = tgt.get("target", "unknown")
                length_histogram = tgt.get("length_histogram", [])

                # Skip targets with no hits
                if not length_histogram:
                    continue

                for length, count in length_histogram:
                    rows.append({
                        "sample": sample_name,
                        "target": target_name,
                        "length": length,
                        "count": count
                    })
        else:
            # Use specific target's histogram
            for tgt in sample.get("targets", []):
                target_name = tgt.get("target", "unknown")
                if target_name == target or target in target_name:
                    length_histogram = tgt.get("length_histogram", [])

                    for length, count in length_histogram:
                        rows.append({
                            "sample": sample_name,
                            "target": target_name,
                            "length": length,
                            "count": count
                        })

    return pd.DataFrame(rows)


def list_targets(input_file):
    """List available targets from the JSON file."""
    with open(input_file, 'r') as f:
        data = json.load(f)

    targets = set()
    for sample in data.get("samples", []):
        for tgt in sample.get("targets", []):
            target_name = tgt.get("target", "unknown")
            reads = tgt.get("reads_with_hits", 0)
            targets.add((target_name, reads))

    return sorted(targets, key=lambda x: (-x[1], x[0]))  # Sort by reads desc, then name


def main():
    parser = argparse.ArgumentParser(
        description="Plot read length distributions from Grate len output"
    )
    parser.add_argument("input_file", help="Input JSON file from 'grate len'")
    parser.add_argument("-o", "--output", help="Output PNG filename (default: <input_prefix>-lenhist.png)")
    parser.add_argument("-f", "--force", action="store_true", help="Overwrite existing output files")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug output")
    parser.add_argument("-t", "--title", default=None,
                        help="Plot title (default: auto-generated)")
    parser.add_argument("--log-scale", action="store_true",
                        help="Use log scale for y-axis")
    parser.add_argument("--max-len", type=int, default=None,
                        help="Maximum x-axis limit")
    parser.add_argument("--facet", action="store_true",
                        help="Facet by sample instead of overlaying")
    parser.add_argument("-b", "--bandwidth", type=float, default=100.0,
                        help="Bandwidth for kernel density estimation smoothing (default: 100.0)")
    parser.add_argument("--target", default="total",
                        help="Which histogram to plot: 'total' (default), 'all', or target name")
    parser.add_argument("--list-targets", action="store_true",
                        help="List available targets and exit")
    parser.add_argument("--by-target", action="store_true",
                        help="When using --target=all, facet by target instead of overlaying")

    args = parser.parse_args()

    # Check if input file exists first
    if not os.path.exists(args.input_file):
        print(f"ERROR: File {args.input_file} does not exist")
        sys.exit(1)

    # Handle --list-targets
    if args.list_targets:
        targets = list_targets(args.input_file)
        print("Available targets:")
        for name, reads in targets:
            print(f"  {name}: {reads} reads")
        sys.exit(0)

    # Auto-generate output filename from input prefix if not specified
    input_filename = os.path.basename(args.input_file)
    input_prefix = os.path.splitext(input_filename)[0]

    if args.output is None:
        suffix = f"-{args.target}" if args.target != "total" else ""
        args.output = f"{input_prefix}{suffix}-lenhist.png"

    # Check if output file exists and require --force to overwrite
    if os.path.exists(args.output) and not args.force:
        print(f"ERROR: Output file already exists: {args.output}")
        print("Use --force to overwrite existing files")
        sys.exit(1)

    try:
        if os.path.getsize(args.input_file) == 0:
            print(f"ERROR: File {args.input_file} is empty")
            sys.exit(1)

        df = load_data(args.input_file, target=args.target)

        if args.debug:
            print(f"\n{args.input_file}:")
            print(f"  Columns: {list(df.columns)}")
            print(f"  Shape: {df.shape}")
            print("  First few rows:")
            print(df.head())

        if df.empty:
            print("ERROR: Input file has no data (or no matching targets)")
            sys.exit(1)

        # Ensure required columns exist
        required_cols = {"sample", "length", "count"}
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            print(f"ERROR: Data missing required columns: {', '.join(missing)}")
            sys.exit(1)

        # Expand data for KDE: repeat each length by its count
        expanded_rows = []
        for _, row in df.iterrows():
            for _ in range(int(row["count"])):
                expanded_rows.append({
                    "sample": row["sample"],
                    "target": row.get("target", "Total"),
                    "length": row["length"]
                })

        kde_df = pd.DataFrame(expanded_rows)

        # Determine groupby columns based on mode
        if args.target == "all":
            groupby_cols = ["sample", "target"]
            color_field = "target:N"
            color_title = "Target"
        else:
            groupby_cols = ["sample"]
            color_field = "sample:N"
            color_title = "Sample"

        # Add total count per group for scaling density to frequency
        kde_df['total_count'] = kde_df.groupby(groupby_cols)['length'].transform('count')

        if args.debug:
            print(f"\nExpanded data for KDE:")
            print(f"  Shape: {kde_df.shape}")
            print(f"  Total reads: {len(kde_df)}")
            print(f"  Length range: {kde_df['length'].min()} to {kde_df['length'].max()}")
            print(f"  Samples: {kde_df['sample'].unique().tolist()}")
            if "target" in kde_df.columns:
                print(f"  Targets: {kde_df['target'].unique().tolist()}")

        # Ordering: preserve data order (order of first appearance in data file)
        sample_order = df["sample"].dropna().astype(str).unique().tolist()
        target_order = df["target"].dropna().astype(str).unique().tolist() if "target" in df.columns else []

        # Auto-generate title
        if args.title is None:
            if args.target == "total":
                args.title = "Read length distribution (all targets)"
            elif args.target == "all":
                args.title = "Read length distribution by target"
            else:
                args.title = f"Read length distribution: {args.target}"

        # --- Plot ---
        alt.data_transformers.enable("json")

        # Determine y-axis scale
        y_scale = alt.Scale(type="log") if args.log_scale else alt.Scale(type="linear")

        # Determine x-axis scale
        x_scale_kwargs = {"zero": False}
        if args.max_len is not None:
            x_scale_kwargs["domain"] = [0, args.max_len]

        # Build groupby for transform_density
        density_groupby = groupby_cols + ["total_count"]

        # Base chart with KDE
        base = (
            alt.Chart(kde_df)
            .transform_density(
                density="length",
                groupby=density_groupby,
                bandwidth=args.bandwidth,
                as_=["length", "density"]
            )
            .transform_calculate(
                frequency="datum.density * datum.total_count"
            )
            .mark_line(strokeWidth=2, interpolate="monotone", clip=True)
            .encode(
                x=alt.X("length:Q",
                       title="Read length (bp)",
                       scale=alt.Scale(**x_scale_kwargs)),
                y=alt.Y("frequency:Q",
                       title="Frequency",
                       scale=y_scale),
                color=alt.Color(color_field,
                              title=color_title,
                              sort=target_order if args.target == "all" else sample_order,
                              scale=alt.Scale(scheme="category20")),
                tooltip=[
                    alt.Tooltip("sample:N", title="Sample"),
                    alt.Tooltip("target:N", title="Target") if "target" in kde_df.columns else None,
                    alt.Tooltip("length:Q", title="Length (bp)", format=".1f"),
                    alt.Tooltip("frequency:Q", title="Frequency", format=",.1f"),
                ],
            )
        )

        # Filter out None tooltips
        base = base.encode(tooltip=[t for t in base.encoding.tooltip if t is not None])

        if args.facet:
            # Faceted by sample
            chart = (
                base.properties(
                    width=700,
                    height=200
                )
                .facet(
                    row=alt.Row("sample:N",
                               title=None,
                               sort=sample_order,
                               header=alt.Header(labelAlign="left", labelAnchor="start", labelFontSize=12))
                )
                .resolve_scale(y="independent")
                .properties(title=args.title)
                .configure_legend(titleFontSize=12, labelFontSize=11)
                .configure_axis(labelFontSize=11, titleFontSize=12)
                .configure_title(fontSize=14)
            )
        elif args.by_target and args.target == "all":
            # Faceted by target
            chart = (
                base.properties(
                    width=700,
                    height=150
                )
                .facet(
                    row=alt.Row("target:N",
                               title=None,
                               sort=target_order,
                               header=alt.Header(labelAlign="left", labelAnchor="start", labelFontSize=11))
                )
                .resolve_scale(y="independent")
                .properties(title=args.title)
                .configure_legend(titleFontSize=12, labelFontSize=11)
                .configure_axis(labelFontSize=11, titleFontSize=12)
                .configure_title(fontSize=14)
            )
        else:
            # Overlay all samples on one plot
            chart = (
                base.properties(
                    width=800,
                    height=500,
                    title=args.title
                )
                .configure_legend(titleFontSize=12, labelFontSize=11)
                .configure_axis(labelFontSize=11, titleFontSize=12)
                .configure_title(fontSize=14)
            )

        chart.save(args.output, scale_factor=2.0)
        print(f"Plot saved to: {args.output}")

    except Exception as e:
        print(f"Error creating plot: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
