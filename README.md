# nanorust

Fast Rust reimplementation of the nanoT BrdU parsing pipeline (steps 01–04 of
`parsing_DoradoRemora_v18_Br.r` + step 02 `bamCoverage`):

```
mod_calls_<exp>.bam  ──►  <prefix>_nanoT_alldata.rds   (same tibble as step 04 `alldata`)
                     └─►  <prefix>.bw                  (same values as deeptools 3.5.4 bamCoverage)
```

One pass over the unsorted dorado BAM. No splitting, no sorting, no samtools/R/deeptools.

## Install

**Prebuilt binaries** (no dependencies; Linux binaries are static and run on any distribution):
download the archive for your platform from the
[Releases page](https://github.com/touala/nanorust/releases), then

```bash
tar xzf nanorust-*.tar.gz && ./nanorust --help
```

**From source** (needs a Rust toolchain from <https://rustup.rs> and a C compiler, Linux or macOS):

```bash
cargo install --git https://github.com/touala/nanorust
```

## Run

```bash
nanorust run -b mod_calls_FP27_ES1_tel1rif1rif2_RefBT1mono_modl610FT11.bam \
             -o FP27_ES1_tel1rif1rif2_RefBT1mono_modl610FT11 -t 16
```

| option | default | meaning |
|---|---|---|
| `-o, --out-prefix` | | writes `<prefix>_nanoT_alldata.rds` and `<prefix>.bw` |
| `--rds`, `--bw` | | explicit output paths (either can be omitted) |
| `-t, --threads` | all cores | **total** CPU cores used (see below) |
| `--no-supplementary` | off | drop supplementary mappings (R keeps same-chrom/strand ones) |
| `--max-dist` | 15000 | `supp_filter(max_dist=)` |
| `--min-len` | 1 | `extract.local.signal(min_len=)` |
| `--bin-size` | 1000 | signalbin bin size |
| `--cov-bin-size` | 50 | bamCoverage `--binSize` |
| `--exclude-prefix` | chrM | chromosomes excluded from the signal (step 01 `^chrM`) |
| `--preset` | brdu | `brdu` (BrdU) or `bredu` (BrdU + EdU), see below |
| `--mod` | | override the preset with one or two modifications: `B`, `E`, or MM specs like `T+B,T+E` (same base) |
| `--source` | auto | supplementary rules: `auto` (from the MM tags), `dorado`, `dnascent` |
| `--binarise THR` | off | R's `binarise`/`bin_thr` for every modification: `signal<X>` uses `prob < THR ? 0 : 1` |
| `--float-mode` | x87 | `x87` = R on x86-64 Linux (80-bit long double), `f64` = R on arm64 macOS |

### Modifications and BAM sources

The producer of the BAM is recognised from its `MM` tags, and the matching R parser is reproduced:

| BAM | MM entries | supplementary rules | `--preset brdu` | `--preset bredu` |
|---|---|---|---|---|
| dorado / Remora | `T+B.` `T+E.` | flag 0/16, or 2048/2064 with first `SA` entry on the same chrom/strand, then `supp_filter` | `parsing_DoradoRemora_v18_Br.r` | `parsing_DoradoRemora_v18_BrE.r` |
| DNAscent 4 | `N+b?` `N+e?` | all mapped, non-secondary records; no `supp_filter` | — | `parsing_DS421_v2.r` |

Output columns: with one modification X, `signalbin(positions, signalX)`, `med_signal`, `med_signalbin`
(v18 layout); with two modifications X and Y, `signalbin(positions, signalX, signalY)`,
`med_signalXbin`, `med_signalYbin` (BrE / DS421 layout). `B`/`b` and `E`/`e` are custom codes of these
models (not in the SAMtags table), so both spellings map to the `B` and `E` columns.

With two modifications, a position is kept when either one has a value (in these BAMs both always list the
same positions). Otherwise the missing value is `NA`, an empty bin mean is `NaN` and an empty median is
`NA`, as in R, and a warning reports how many mappings were affected. Both modifications must be on the
same base.

### CPU budget (`-t`)

`-t N` is the total number of cores the run uses, so it can be set to the CPUs allocated by
HTCondor (`request_cpus`) or SLURM (`--cpus-per-task`):

* `-t 1`: everything runs in one thread.
* `-t 2`: one thread reads and decompresses, one processes records.
* `-t N` (N ≥ 3): decompression and processing share N−1 compute slots; the remaining core is
  left to the light file-reading and record-splitting threads. Output (RDS compression, then
  bigWig) runs after processing, within the same N.

Measured average busy cores on a 3.8 GB BAM: 0.99 (`-t 1`), 1.5 (`-t 2`), 2.3 (`-t 3`),
3.5 (`-t 4`), 5.2 (`-t 6`). Results are identical for any `-t`.

`-b -` reads the BAM from stdin, e.g. streaming from a server without a local copy:

```bash
ssh server cat /path/mod_calls.bam | nanorust run -b - -o PREFIX
```

`nanorust bw-compare a.bw b.bw` compares two bigWigs base by base.

## What is replicated

* **Filters** – mapped, non-secondary, not `chrM*`, has `MM`/`ML`; flag 0/16, or 2048/2064 when the
  first `SA` entry is on the same chrom and strand; `end - start > min_len` with `rlen = M + D`.
* **Mod tags** – the selected calls over the target bases of the read-oriented sequence (complement on
  minus reads; listed positions for DNAscent's `N`), `.` → unreported = 0, `?` → dropped, `ML/255`;
  combined entries such as `C+mh` are supported.
* **CIGAR mapping** – `parseCigar` semantics (M maps, I/S advance query, D/N advance reference,
  minus-strand flip using `max(read_pos)`), positions kept in `[start, end]`.
* **Binning** – `floor((pos-1)/1000)*1000+1`, mean per bin; `med_signal`, `med_signalbin` with R's
  `mean`/`median` algorithms.
* **supp_filter** – per read: total-or-null read-position overlap with the first mapping and
  distance < `max_dist`.
* **RDS** – R serialization v3, gzip; tibble with the same columns, types, factor levels
  (all BAM references, strand `+ - *`) and nested `signalbin` tibbles.
* **Coverage** – deeptools 3.5.4 `bamCoverage` defaults on `samtools view -F 260`: per 50 bp bin,
  number of reads with an aligned block (pysam `get_blocks`) in the bin.

## Known, intended differences to the R output

* **Row order**: rows are sorted by `(chrom, read_id, flag, start)` over the whole file, ties in BAM
  order (R's stable `arrange()`); R orders within each 10k-read chunk, so the order matches R exactly
  when R processed the BAM as a single chunk. `supp_filter` also runs per whole read instead of per chunk.
* **Floating point**: R's `mean()` accumulates in `long double`, which is 80-bit x87 on x86-64
  Linux. nanorust emulates it in software (default `--float-mode x87`), so values are bit-identical
  to R on Linux on any machine; `--float-mode f64` reproduces R on Apple Silicon instead.
* **Downstream row-order effect**: step 04's `mean`/`var` sum `signalB` in `alldata` row order, so
  with nanorust's order `mean_br_bin`/`varbin` can differ from the R pipeline in the last bit
  (≤ 4e-16); `nbin` and the median filter are identical. With rows in R's order, step 04 is identical.
* The bigWig has the same values/intervals but is not byte-identical (different writer library).

## Validation

```bash
scripts/validate.sh sample/chrI.bam ref_nanoT_alldata.rds ref.bw chrI
```

* Whole genome, FP27_ES1_tel1rif1rif2 (1,296,021 records, 11 GB): the same 313,924 mappings as R;
  coverage bigWig identical at every base on all 17 chromosomes (same intervals).
* chrI, chrVI, chrXII (105,307 mappings) with `--float-mode x87`: `alldata` is `identical()` to the
  Linux R output after ordering, every value bit for bit.
* BrdU + EdU (`--preset bredu`) on human CHM13 test sets, dorado and DNAscent, 1,000 and 10,000 reads
  (18,781 mappings): the whole `alldata` is `identical()` to the R parser output (minus its step-03
  `signal` column), row order included.
* Other modifications, combined MM codes and edge cases: end-to-end tests on rewritten BAMs.

## License

MIT
