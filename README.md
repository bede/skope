[![Crates.io version](https://img.shields.io/crates/v/skope?style=flat-square)](https://crates.io/crates/skope)
[![Conda version](https://img.shields.io/conda/v/bioconda/skope?style=flat-square&label=bioconda&color=blue)](https://anaconda.org/bioconda/skope)

# Skope

Accelerated streaming containment and abundance estimation using syncmers. Like Mash Screen or Sylph, only much faster for small queries (up to 10 gigabases/second), despite more dense and even *k*-mer sampling. Rapidly estimates coverage at chosen depth thresholds, providing a viable alternative to mapping-based coverage estimation in many contexts. Optimised for screening large datasets for small query sequences (<1 gigabase) such as genes and genomes of viruses and bacteria, scaling to larger query genomes given sufficient memory. No prior sketching is required, and memory use is bounded by the size of query (needle) sequences rather than by the size of subject (haystack) sequences.

## Install & update

```bash
# Bioconda
conda install -c bioconda deacon

# Latest stable
RUSTFLAGS="-C target-cpu=native" cargo install skope

# Latest from git
RUSTFLAGS="-C target-cpu=native" cargo install --git https://github.com/bede/skope
```

## Usage

```bash
# Calculate containment of target sequences in reads
skope query refs.fa reads.fastq.gz

# Treat each record in a multi-record fastx as a separate target
skope query -i refs.fa reads.fastq.gz

# Bypass syncmer selection and take every canonical k-mer
# (needed when targets are bare k-mers: a 31 bp record survives selection only ~1 time in 23)
skope query --all-kmers -i kmers31.fa reads.fastq.gz

# Calculate target containment in multiple samples
skope query refs.fa reads1.fastq.gz reads2/ reads3.fa.zst…

# Calculate containment at depth>=100 using only discriminatory k-mers
skope query -i -a 100 --discriminatory refs.fa reads1.fastq.gz reads2/ reads3.fa.zst…

# Plot query results (containment bar chart or scatter)
skope query -i refs.fa s1.fq.gz s2.fq.gz > query.tsv
uv run plot/query.py query.tsv --mode bar -o query-bar.png
uv run plot/query.py query.tsv --mode decay -o query-decay.png
uv run plot/query.py query.tsv --mode scatter -o query-scatter.png
uv run plot/query.py query.tsv --mode scatter -o query-scatter.html  # Interactive

# Stdin
zstdcat reads3.fq.zst | skope query refs.fa -

# Mask off-target background k-mers in background.fa
skope query refs.fa -b background.fa reads.fq.gz
# …or baked once into a reusable query index (.sk)
skope index build-query refs.fa -b background.fa -o refs.sk
skope query refs.sk reads.fq.gz

# The build commands take targets on stdin too (they have no samples competing for it)
zstdcat refs.fa.zst | skope index build-query - -o refs.sk
```

Run the plotting scripts with [uv](https://docs.astral.sh/uv/) to automatically handle dependencies.

**Plotting scripts:**
- `plot/query.py` - Containment bar charts and scatter plots from `query` TSV output

![Example containment plot](data/multi.png)

### CLI Reference

**Main commands**
```
  query     Estimate k-mer containment & abundance in fastx file(s) or directories
  index     Build and manage query and classification indexes (alpha)
```

**Query**

```bash
$ skope query -h
Estimate target containment & abundance in fastx file(s) or directories thereof using open syncmers or all k-mers

Usage: skope query [OPTIONS] <TARGETS> <SAMPLES>...

Arguments:
  <TARGETS>     Path to fastx file (single target unless -i), directory of fastx files/subdirs (one target per child file/subdir) or query index (.sk)
  <SAMPLES>...  Path(s) to fastx files/dirs (- for stdin). Each file/dir is treated as a separate sample

Options:
  -k, --kmer <K>                        K-mer length (1-61, default 31), read from a prebuilt index
  -s, --smer <S>                        S-mer length (odd, s < k, default 9), read from a prebuilt index
  -i, --individual                      Treat each fastx record as separate target (default: merge records into one target)
  -c, --confidence                      Report confidence intervals, ANI estimates, and patchiness columns
  -d, --discriminatory                  Consider only k-mers unique to each target
      --all-kmers                       Evaluate all k-mers, bypassing syncmer selection
      --complexity <FLOAT>              FracMinHash fraction of targeters below this kdust complexity [0, 1] (0 = retain all) [default: 0]
  -f, --fraction <FLOAT>                Fraction of target k-mers to keep [0, 1] [default: 1]
  -a, --abundance-thresholds <INT,...>  Comma-separated additional abundance thresholds for containment estimation [default: 10]
  -b, --background <BACKGROUND>         Path to fastx file(s) whose k-mers we wish to drop from our targets
  -l, --limit <BASES>                   Terminate processing after approximately this many bases (e.g. 50M, 10G)
  -t, --threads <THREADS>               Number of execution threads (0 = auto) [default: 8]
  -o, --output <OUTPUT>                 Path to output file (- for stdout) [default: -]
  -n, --names <NAME,...>                Comma-separated sample names (default is file/dir name without extension)
      --sort <SORT>                     Sort results [default: containment] [possible values: containment, target, input]
      --dump-kmers <FILE>               Dump selected target k-mers to TSV file (target, position, kmer)
      --no-total                        Suppress TOTAL summary rows in output
  -q, --quiet                           Suppress progress reporting
  -h, --help                            Print help
```

## Confidence

Passing `--confidence` (`-c`) to `skope query` adds output columns for confidence, estimated ANI, and coverage 'patchiness'.

```bash
skope query --confidence refs.fa reads.fq
```

- `containmentX_ci`: a 95% Wilson score confidence interval for each containment estimate reflecting uncertainty in the proportion `hits/target_kmers` at each abundance threshold. This is written as `{lower}-{upper}`. Intervals are narrower for long target sequences and wider for short targets and/or low containment.

- `patchiness`: a Wald–Wolfowitz runs test for clustering of `containment1` hits along the target sequence. Written as `{z}|{p}`, positive `z` means that selected *k*-mer distribution is more more patchy than expected for the same hit count, and `p` is the one-sided normal-approximation p-value. Displayed only for `z > 0` and `p <= 0.05` (otherwise `-`), which can also mean the test was skipped for too few eligible k-mers, hits, or misses.

- `ani_est`: a containment ANI estimate based on `containment1`. Skope transforms containment with `containment^(1/k)`, and when the low-coverage depth histogram is sufficient it first applies the Sylph-like Poisson sampling adjustment `containment1 / Pr(Pois(λ) >= 1)`. If the adjustment cannot be estimated, the unadjusted containment ANI is reported instead.

  `ani_est` is shown as `-` (suppressed) for any target with fewer than 50 k-mers, with no contained k-mers, or whose estimate falls below 0.90, since these yield too little signal for a meaningful estimate.

  Calculating `ani_est` is skipped under the following conditions:

  - Median nonzero depth > 2

  - Fewer than 50 target k-mers

  - Fewer than 25 hitting k-mers

  - Fewer than 3 k-mers at either the modal nonzero depth or the depth above it

## Dumping k-mers

Passing `--dump-kmers <path>` to `skope query` writes the selected target k-mers to a TSV file with columns `target`, `position`, and `kmer`. The dump reflects whatever selection is in effect, so `--discriminatory` is respected.

One row is emitted per k-mer occurrence, so the row count can exceed `target_kmers` (which counts distinct k-mers) when a *k*-mer recurs within a target. The `kmer` column is the canonical *k*-mer, not necessarily the forward-strand sequence at that position.
