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


def main():
    parser = argparse.ArgumentParser(description="Plot Skope TSV file (single TSV; one or many samples)")
    # Input/output
    parser.add_argument("input_tsv", help="Input TSV file containing one or many samples (must include a sample column)")
    parser.add_argument("-o", "--output", help="Output plot filename (default: <input_prefix>.png)")
    parser.add_argument("-f", "--force", action="store_true", help="Overwrite existing output files")
    parser.add_argument("--format", choices=["png", "html"], default="png",
                        help="Output format (default: %(default)s)")
    # Data selection
    parser.add_argument("-s", "--samples", help="Only plot these samples (comma-separated)")
    parser.add_argument("--sample-column", default="sample",
                        help="Name of the column that identifies samples (default: 'sample')")
    parser.add_argument("--totals", action="store_true",
                        help="Plot only TOTAL rows (aggregate stats per sample)")
    parser.add_argument("--totals-label",
                        help="Custom label to replace 'TOTAL' when plotting totals")
    parser.add_argument("-a", "--abundance-threshold", type=int, default=1,
                        help="Abundance threshold for containment column to plot (default: %(default)s for containment1)")
    # Plot mode
    parser.add_argument("-m", "--mode", choices=["bar", "scatter"], default="bar",
                        help="Plot mode: 'bar' for bar chart, 'scatter' for containment vs abundance (default: %(default)s)")
    parser.add_argument("--log-y", action="store_true",
                        help="Use logarithmic scale for y-axis (scatter mode only)")
    # Appearance
    parser.add_argument("-t", "--title", default="Containment analysis",
                        help="Plot title (default: %(default)s)")
    parser.add_argument("--short-names", action="store_true",
                        help="Remove accession prefix (before first space) from target names")
    parser.add_argument("-c", "--colours", default="category20",
                        help="Altair colour scheme (default: %(default)s)")
    parser.add_argument("--no-depth", action="store_true",
                        help="Disable depth labels on the plot")
    # Debug
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug output")

    args = parser.parse_args()

    # Auto-generate output filenames from input prefix if not specified
    # Use only the filename (not path) and output to current directory
    input_filename = os.path.basename(args.input_tsv)
    input_prefix = os.path.splitext(input_filename)[0]
    if args.output is None:
        args.output = f"{input_prefix}.{args.format}"
    else:
        # Infer format from output extension
        ext = os.path.splitext(args.output)[1].lower()
        if ext == ".html":
            args.format = "html"
        elif ext == ".png":
            args.format = "png"

    # Check if output files exist and require --force to overwrite
    if os.path.exists(args.output) and not args.force:
        print(f"ERROR: Output file already exists: {args.output}")
        print("Use --force to overwrite existing file")
        sys.exit(1)

    try:
        # --- Load & validate ---
        if not os.path.exists(args.input_tsv):
            print(f"ERROR: File {args.input_tsv} does not exist")
            sys.exit(1)
        if os.path.getsize(args.input_tsv) == 0:
            print(f"ERROR: File {args.input_tsv} is empty")
            sys.exit(1)

        df = pd.read_csv(args.input_tsv, sep='\t')

        if args.debug:
            print(f"\n{args.input_tsv}:")
            print(f"  Columns: {list(df.columns)}")
            print(f"  Shape: {df.shape}")
            print("  First few rows:")
            print(df.head())

        if df.empty:
            print("ERROR: Input TSV has no data")
            sys.exit(1)

        # Construct containment column name from abundance threshold
        containment_col = f"containment{args.abundance_threshold}"
        hits_col = f"containment{args.abundance_threshold}_hits"

        # Ensure required columns exist
        required_cols = {"target", containment_col, hits_col, "median_nz_abundance"}
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            print(f"ERROR: Input TSV missing required columns: {', '.join(missing)}")
            sys.exit(1)

        sample_col = args.sample_column
        if sample_col not in df.columns:
            print(f"ERROR: Input TSV missing sample column '{sample_col}'. "
                  f"Use --sample-column to specify the correct column.")
            sys.exit(1)

        # Coerce numeric columns in case they're strings
        for c in (containment_col, hits_col, "median_nz_abundance", "length_bp", "contained_minimizers"):
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")

        # --- Prep for plotting ---
        plot_df = df.copy()

        # Filter TOTAL rows based on --totals flag
        if "target" in plot_df.columns:
            if args.totals:
                # Keep only TOTAL rows
                plot_df = plot_df[plot_df["target"] == "TOTAL"]
            else:
                # Drop TOTAL rows (default behavior)
                plot_df = plot_df[plot_df["target"] != "TOTAL"]

        if plot_df.empty:
            print("ERROR: No data to plot after filtering")
            sys.exit(1)

        # Filter to specific samples if requested
        if args.samples:
            samples_list = [s.strip() for s in args.samples.split(",")]
            plot_df = plot_df[plot_df[sample_col].isin(samples_list)]
            if plot_df.empty:
                print(f"ERROR: No matching samples found. Available: {df[sample_col].unique().tolist()}")
                sys.exit(1)

        # Display name tweaks
        if args.short_names:
            plot_df["display_name"] = plot_df["target"].apply(
                lambda x: str(x).split(" ", 1)[1] if isinstance(x, str) and " " in x else x
            )
        else:
            plot_df["display_name"] = plot_df["target"]

        # Replace TOTAL label if custom label provided
        if args.totals_label:
            plot_df["display_name"] = plot_df["display_name"].replace("TOTAL", args.totals_label)

        # Ordering: preserve data order (order of first appearance in CSV)
        # Use pandas .unique() which preserves order of appearance
        sample_order = plot_df[sample_col].dropna().astype(str).unique().tolist()
        target_order = plot_df["display_name"].dropna().astype(str).unique().tolist()

        # Sort dataframe by target order for consistent bar grouping
        plot_df["target_idx"] = plot_df["display_name"].apply(lambda x: target_order.index(x) if x in target_order else -1)
        plot_df = plot_df.sort_values(["target_idx", sample_col])
        plot_df = plot_df.drop("target_idx", axis=1)

        # Text labels for median abundance
        plot_df["depth_label"] = plot_df["median_nz_abundance"].apply(
            lambda x: f"med(abund): {x:.0f}" if pd.notna(x) and float(x) > 0 else ""
        )

        # --- Plot ---
        alt.data_transformers.enable("json")

        if args.mode == "scatter":
            # Scatter plot: containment vs abundance
            y_scale = alt.Scale(type="log") if args.log_y else alt.Scale()
            selection = alt.selection_point(fields=[sample_col], bind="legend")
            chart = (
                alt.Chart(plot_df)
                .mark_circle(size=60)
                .encode(
                    x=alt.X(f"{containment_col}:Q", title=f"Containment at depth ≥ {args.abundance_threshold}", scale=alt.Scale(domain=[0, 1])),
                    y=alt.Y("median_nz_abundance:Q", title="Median non-zero abundance", scale=y_scale),
                    color=alt.Color(
                        f"{sample_col}:N",
                        sort=sample_order,
                        title="",
                        scale=alt.Scale(scheme=args.colours),
                    ),
                    opacity=alt.condition(selection, alt.value(1), alt.value(0.1)),
                    tooltip=[
                        alt.Tooltip("target:N", title="Target"),
                        alt.Tooltip(f"{sample_col}:N", title="Sample"),
                        alt.Tooltip("length_bp:Q", title="Length", format=",.0f"),
                        alt.Tooltip(f"{containment_col}:Q", title="Containment"),
                        alt.Tooltip("median_nz_abundance:Q", title="Median abundance"),
                    ],
                )
                .add_params(selection)
                .properties(title=args.title, width=450, height=400)
                .configure_legend(titleFontSize=14, labelFontSize=12, symbolSize=150)
                .configure_axis(labelFontSize=12, titleFontSize=14)
                .configure_title(fontSize=14)
            )
        else:
            # Bar chart (default)
            selection = alt.selection_point(fields=[sample_col], bind="legend")
            bars = (
                alt.Chart(plot_df)
                .mark_bar(size=6)
                .encode(
                    y=alt.Y("display_name:N", title="", sort=target_order),
                    x=alt.X(f"{containment_col}:Q", title=f"Containment at depth ≥ {args.abundance_threshold}", scale=alt.Scale(domain=[0, 1])),
                    color=alt.Color(
                        f"{sample_col}:N",
                        sort=sample_order,
                        title="",
                        scale=alt.Scale(scheme=args.colours),
                    ),
                    yOffset=alt.YOffset(f"{sample_col}:N", sort=sample_order),
                    opacity=alt.condition(selection, alt.value(1), alt.value(0.1)),
                    tooltip=[
                        alt.Tooltip("target:N", title="Target"),
                        alt.Tooltip(f"{sample_col}:N", title="Sample"),
                        alt.Tooltip("length_bp:Q", title="Length", format=",.0f"),
                        alt.Tooltip(f"{containment_col}:Q", title="Containment"),
                        alt.Tooltip("median_nz_abundance:Q", title="Median abundance"),
                    ],
                )
                .add_params(selection)
            )

            # Conditionally add depth labels
            if not args.no_depth:
                text_labels = (
                    alt.Chart(plot_df)
                    .mark_text(align="left", baseline="middle", dx=5, fontSize=7, color="black")
                    .encode(
                        y=alt.Y("display_name:N", sort=target_order),
                        yOffset=alt.YOffset(f"{sample_col}:N", sort=sample_order),
                        x=alt.value(0),
                        text="depth_label:N",
                    )
                )
                chart = bars + text_labels
            else:
                chart = bars

            chart = (
                chart
                .properties(title=args.title, width=450, height=alt.Step(7))
                .configure_legend(titleFontSize=14, labelFontSize=12, symbolSize=150)
                .configure_axis(labelFontSize=12, titleFontSize=14)
                .configure_title(fontSize=14)
                .resolve_scale(y="shared")
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
