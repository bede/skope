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
    parser = argparse.ArgumentParser(
        description="Plot read length distributions from Grate len output"
    )
    parser.add_argument("input_files", nargs="+", help="Input TSV file(s) from 'grate len'")
    parser.add_argument("-o", "--output", help="Output PNG filename (default: <input_prefix>-lenhist.png)")
    parser.add_argument("-f", "--force", action="store_true", help="Overwrite existing output files")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug output")
    parser.add_argument("-t", "--title", default="Read length distribution",
                        help="Plot title (default: %(default)s)")
    parser.add_argument("--log-scale", action="store_true",
                        help="Use log scale for y-axis")
    parser.add_argument("--max-len", type=int, default=None,
                        help="Maximum x-axis limit")
    parser.add_argument("--min-len", type=int, default=0,
                        help="Minimum x-axis limit (default: 0)")
    parser.add_argument("--facet", action="store_true",
                        help="Facet by sample instead of overlaying")
    parser.add_argument("-b", "--bin-step", type=int, default=100,
                        help="Bin step size in bp (default: 100)")
    parser.add_argument("-c", "--colours", "--colors", default="category20",
                        help="Altair/Vega-Lite color scheme (default: %(default)s)")

    args = parser.parse_args()

    input_filename = os.path.basename(args.input_files[0])
    input_prefix = os.path.splitext(input_filename)[0]

    if args.output is None:
        args.output = f"{input_prefix}-lenhist.png"

    if os.path.exists(args.output) and not args.force:
        print(f"ERROR: Output file already exists: {args.output}")
        print("Use --force to overwrite existing files")
        sys.exit(1)

    try:
        dfs = []
        for input_file in args.input_files:
            if not os.path.exists(input_file):
                print(f"ERROR: File {input_file} does not exist")
                sys.exit(1)

            if os.path.getsize(input_file) == 0:
                print(f"ERROR: File {input_file} is empty")
                sys.exit(1)

            file_df = pd.read_csv(input_file, sep='\t')

            if args.debug:
                print(f"\n{input_file}:")
                print(f"  Columns: {list(file_df.columns)}")
                print(f"  Shape: {file_df.shape}")
                print("  First few rows:")
                print(file_df.head())

            dfs.append(file_df)

        df = pd.concat(dfs, ignore_index=True)

        if df.empty:
            print("ERROR: Input files have no data")
            sys.exit(1)

        required_cols = {"sample", "length", "count"}
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            print(f"ERROR: Data missing required columns: {', '.join(missing)}")
            sys.exit(1)

        if args.debug:
            print(f"\nData summary:")
            print(f"  Shape: {df.shape}")
            print(f"  Total reads: {df['count'].sum()}")
            print(f"  Length range: {df['length'].min()} to {df['length'].max()}")
            print(f"  Samples: {df['sample'].unique().tolist()}")

        sample_order = df["sample"].dropna().astype(str).unique().tolist()

        alt.data_transformers.enable("json")

        y_scale = alt.Scale(type="log") if args.log_scale else alt.Scale(type="linear")

        min_len = args.min_len
        max_len = args.max_len if args.max_len else int(df["length"].max())
        bin_extent = [min_len, max_len]

        base = (
            alt.Chart(df)
            .mark_line(interpolate="step-after", strokeWidth=2)
            .encode(
                x=alt.X(
                    "length:Q",
                    bin=alt.Bin(step=args.bin_step, extent=bin_extent),
                    title="Read length (bp)",
                ),
                y=alt.Y(
                    "sum(count):Q",
                    stack=None,
                    title="Read count",
                    scale=y_scale,
                ),
                color=alt.Color(
                    "sample:N",
                    title="Sample",
                    sort=sample_order,
                    scale=alt.Scale(scheme=args.colours),
                ),
            )
        )

        if args.facet:
            chart = (
                base.properties(width=700, height=200)
                .facet(
                    row=alt.Row(
                        "sample:N",
                        title=None,
                        sort=sample_order,
                        header=alt.Header(labelAlign="left", labelAnchor="start", labelFontSize=12),
                    )
                )
                .resolve_scale(y="independent")
                .properties(title=args.title)
                .configure_legend(titleFontSize=12, labelFontSize=11)
                .configure_axis(labelFontSize=11, titleFontSize=12)
                .configure_title(fontSize=14)
            )
        else:
            chart = (
                base.properties(width=800, height=500, title=args.title)
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
