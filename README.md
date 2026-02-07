# Skope

Fast containment and abundance estimation.

## Install & update

```bash
RUSTFLAGS="-C target-cpu=native" cargo install --git https://github.com/bede/skope
```

## Usage

```bash
# Calculate containment of target sequences in reads
skope query refs.fa reads.fastq.gz

# Calculate target containment in multiple samples
skope query refs.fa reads1.fastq.gz reads2/ reads3.fa.zst…

# Output as table instead of TSV
skope query -f table refs.fa reads.fastq.gz

# Plot query results (containment bar chart or scatter)
skope query refs.fa s1.fq.gz s2.fq.gz > query.tsv
uv run plot/query.py query.tsv --mode bar -o query_bar.png
uv run plot/query.py query.tsv --mode scatter -o query_scatter.png

# Plot sequence length histograms
skope lenhist - s1.fq.gz s2.fq.gz > len.tsv
uv run plot/lenhist.py len.tsv -f -b 500

# Stdin
zstdcat reads3.fq.zst | skope query refs.fa -

# View help for any command
skope query -h
skope lenhist -h
skope classify -h
skope index build -h
```

Run the plotting scripts with [uv](https://docs.astral.sh/uv/) to automatically handle dependencies.

**Plotting scripts:**
- `plot/query.py` - Containment bar charts and scatter plots from `query` TSV output
- `plot/lenhist.py` - Length distribution histograms from `lenhist` TSV output

![Example containment plot](data/multi.png)

### CLI Reference

**Main commands:**
```
  query     Estimate k-mer containment & abundance in fastx file(s) or directories thereof
  classify  Classify sequences into groups based on k-mer membership
  lenhist   Generate length histogram for sequences with k-mer hits to target sequence(s)
  index     Build and manage classification indexes
```

**Query containment:**
```bash
$ skope query -h
Estimate k-mer containment & abundance in fastx file(s) or directories thereof

Usage: skope query [OPTIONS] <TARGETS> <SAMPLES>...

Arguments:
  <TARGETS>     Path to fastx file containing target sequence record(s)
  <SAMPLES>...  Path(s) to fastx files/dirs (- for stdin). Each file/dir is treated as a separate sample

Options:
  -k, --kmer-length <KMER_LENGTH>
          K-mer length (1-61) [default: 31]
  -s, --smer-length <SMER_LENGTH>
          S-mer length used for open syncmer selection (s < k, s must be odd) [default: 9]
  -a, --abundance-thresholds <ABUNDANCE_THRESHOLDS>
          Comma-separated additional abundance thresholds for containment estimation [default: 10]
  -d, --discriminatory
          Consider only syncmers unique to each target
  -t, --threads <THREADS>
          Number of execution threads (0 = auto) [default: 8]
  -l, --limit <LIMIT>
          Terminate processing after approximately this many bases (e.g. 50M, 10G)
  -o, --output <OUTPUT>
          Path to output file (- for stdout) [default: -]
  -f, --format <FORMAT>
          Output format [default: tsv] [possible values: tsv, table]
  -n, --names <SAMPLE_NAMES>
          Comma-separated sample names (default is file/dir name without extension)
  -S, --sort <SORT>
          Sort displayed results: c=containment (descending), t=target, s=sample, o=original [default: c] [possible values: c, t, o, s]
  -q, --quiet
          Suppress progress reporting
      --no-total
          Suppress TOTAL summary rows in output
      --dump-positions <DUMP_POSITIONS>
          Dump open syncmer positions to TSV file (target\tposition)
  -h, --help
          Print help
```

**Length histogram:**
```bash
$ skope lenhist -h
Generate length histogram for sequences with k-mer hits to target sequence(s)

Usage: skope lenhist [OPTIONS] <TARGETS> <SAMPLES>...

Arguments:
  <TARGETS>     Path to fastx file containing target sequence record(s) (- to disable target filtering)
  <SAMPLES>...  Path(s) to fastx files/dirs (- for stdin). Each file/dir is treated as a separate sample

Options:
  -k, --kmer-length <KMER_LENGTH>  k-mer length (1-61) [default: 31]
  -s, --smer-length <SMER_LENGTH>  S-mer length used for open syncmer selection (s < k, s must be odd) [default: 9]
  -t, --threads <THREADS>          Number of execution threads (0 = auto) [default: 8]
  -l, --limit <LIMIT>              Terminate processing after approximately this many bases (e.g. 50M, 10G)
  -o, --output <OUTPUT>            Path to output file (- for stdout) [default: -]
  -n, --names <SAMPLE_NAMES>       Comma-separated sample names (default is file/dir name without extension)
  -q, --quiet                      Suppress progress reporting
  -h, --help                       Print help
```

