#!/usr/bin/env -S uv run
# /// script
# dependencies = [
#   "pandas",
#   "altair",
#   "vl-convert-python",
# ]
# ///

import argparse
import os
import sys

import altair as alt
import pandas as pd


def load_data(input_file):
    """Load length histogram data from CSV file."""
    return pd.read_csv(input_file)


def main():
    parser = argparse.ArgumentParser(
        description="Plot read length distributions from Grate len output"
    )
    parser.add_argument("input_file", help="Input CSV file from 'grate len'")
    parser.add_argument("-o", "--output", help="Output PNG filename (default: <input_prefix>-lenhist.png)")
    parser.add_argument("-f", "--force", action="store_true", help="Overwrite existing output files")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug output")
    parser.add_argument("-t", "--title", default="Read length distribution",
                        help="Plot title (default: %(default)s)")
    parser.add_argument("--log-scale", action="store_true",
                        help="Use log scale for y-axis")
    parser.add_argument("--max-len", type=int, default=None,
                        help="Maximum x-axis limit")
    parser.add_argument("--facet", action="store_true",
                        help="Facet by sample instead of overlaying")
    parser.add_argument("-b", "--bandwidth", type=float, default=100.0,
                        help="Bandwidth for kernel density estimation smoothing (default: 100.0)")

    args = parser.parse_args()

    # Auto-generate output filename from input prefix if not specified
    input_filename = os.path.basename(args.input_file)
    input_prefix = os.path.splitext(input_filename)[0]

    if args.output is None:
        args.output = f"{input_prefix}-lenhist.png"

    # Check if output file exists and require --force to overwrite
    if os.path.exists(args.output) and not args.force:
        print(f"ERROR: Output file already exists: {args.output}")
        print("Use --force to overwrite existing files")
        sys.exit(1)

    try:
        # --- Load & validate ---
        if not os.path.exists(args.input_file):
            print(f"ERROR: File {args.input_file} does not exist")
            sys.exit(1)

        if os.path.getsize(args.input_file) == 0:
            print(f"ERROR: File {args.input_file} is empty")
            sys.exit(1)

        df = load_data(args.input_file)

        if args.debug:
            print(f"\n{args.input_file}:")
            print(f"  Columns: {list(df.columns)}")
            print(f"  Shape: {df.shape}")
            print("  First few rows:")
            print(df.head())

        if df.empty:
            print("ERROR: Input file has no data")
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
                    "length": row["length"]
                })

        kde_df = pd.DataFrame(expanded_rows)

        # Add total count per group for scaling density to frequency
        kde_df['total_count'] = kde_df.groupby(['sample'])['length'].transform('count')

        if args.debug:
            print(f"\nExpanded data for KDE:")
            print(f"  Shape: {kde_df.shape}")
            print(f"  Total reads: {len(kde_df)}")
            print(f"  Length range: {kde_df['length'].min()} to {kde_df['length'].max()}")
            print(f"  Samples: {kde_df['sample'].unique().tolist()}")

        # Ordering: preserve data order (order of first appearance in data file)
        sample_order = df["sample"].dropna().astype(str).unique().tolist()

        # --- Plot ---
        alt.data_transformers.enable("json")

        # Determine y-axis scale
        y_scale = alt.Scale(type="log") if args.log_scale else alt.Scale(type="linear")

        # Determine x-axis scale
        x_scale_kwargs = {"zero": False}
        if args.max_len is not None:
            x_scale_kwargs["domain"] = [0, args.max_len]

        # Base chart with KDE
        base = (
            alt.Chart(kde_df)
            .transform_density(
                density="length",
                groupby=["sample", "total_count"],
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
                color=alt.Color("sample:N",
                              title="Sample",
                              sort=sample_order,
                              scale=alt.Scale(scheme="category20")),
                tooltip=[
                    alt.Tooltip("sample:N", title="Sample"),
                    alt.Tooltip("length:Q", title="Length (bp)", format=".1f"),
                    alt.Tooltip("frequency:Q", title="Frequency", format=",.1f"),
                ],
            )
        )

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
