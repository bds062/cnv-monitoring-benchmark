# cnv-monitoring-benchmark

Code and data tables for the paper **Temporal CNV calling for genomic monitoring** (manuscript in preparation; the citation will be added on publication).

This repository has everything needed to rerun the analysis from the public sequencing reads:
- download scripts;
- one workflow script per CNV caller;
- the scripts that make every figure and table in the paper;
- the small data tables the paper reports, including the ddPCR copy numbers ([`data/ddpcr_tryp77.tsv`](data/ddpcr_tryp77.tsv)).

## What the paper does

A copy-number variant (CNV) is a stretch of DNA present in more or fewer copies than usual. In the corn earworm, *Helicoverpa zea*, an amplification of the trypsin gene *Tryp77* is associated with resistance to Bt toxins. A monitoring program would track how common the amplification is in field populations over time. The paper asks which short-read CNV caller suits that job.

The paper compares five callers. They use two kinds of evidence:
- **Breakpoint callers** find the ends of a CNV from split and discordant reads: Delly2, and smoove (LUMPY + SVTyper).
- **Read-depth callers** compare how many reads fall in each window with other genomes: GATK-gCNV, CNVcaller and CNVkit.

Each caller compares a genome with a set of other genomes, its *background population*. A trio member is never in a background with its own parents or child, except in the pooled-12 comparison.

The callers are scored in three ways:
1. **Mendelian inconsistency in trios.** In four *H. zea* trios and four human 1000 Genomes trios, a CNV called in a child but in neither parent is most likely an error. The share of such calls is the Mendelian inconsistency rate (MIR). It needs no truth set. In human, the calls are also compared with long-read truth sets.
2. **Copy number at *Tryp77*.** The callers' copy numbers in the 12 trio moths are compared with droplet digital PCR (ddPCR).
3. **Monitoring.** The callers are run on published field moths from Louisiana (2002, 2012 and 2017), as a monitoring program would run them.

## Repository layout

```
config.sh       paths and thread count shared by every workflow script
envs/           conda environments, one per tool group (versions pinned)
data/           download.sh, sample sheets, ddPCR copy numbers and other small tables (see data/README.md)
workflow/
  zea/          H. zea: reference, read trimming and alignment, then the five callers
  human/        human: reference, CRAM to BAM, then the five callers
  callers/      one script per caller, shared by both species
  monitoring/   the monitoring analyses on the published Louisiana moths
  downsample/   chromosome 7 of the 12 trio moths at full depth, 10x and 5x, every caller rerun per depth
  runtime/      the fixed-configuration re-timing of the trio benchmark
analysis/       the scripts that make the paper's tables and figures (cnvbench/ holds shared code)
```

Outputs go to `results/` by default:
- Caller outputs: `results/calls/{species}/{caller}/{background}/{sample}/`.
- Tables: `results/tables/`.
- Figures: `results/figures/`, under the file names in the table below.

To keep data or results elsewhere, set `CNVB_DATA` and `CNVB_RESULTS`, either in the environment or in a file `config.local.sh` (not tracked).

## Running the analysis

The steps run in order. Steps 2 and 3 are the long ones: the callers run on 65 genomes (53 moths and 12 people). A compute cluster is recommended.

**1. Install the software.** Create the environments in `envs/`, for example:
```
for f in envs/*.yml; do conda env create -f "$f"; done
```
smoove runs from a container, so `apptainer` is also needed. The figures use Times New Roman if its TTF files are in `~/.fonts`; without them, matplotlib falls back to DejaVu Serif and only the font changes.

**2. Download the data and the remaining tools.** All inputs are public; `data/README.md` lists each one.
```
bash data/download.sh
```
The reads of the four *H. zea* trio children will be released in the NCBI Sequence Read Archive with the paper, and `download.sh` skips them until then. All other reads are already public.

**3. Run the callers.** Run one stage at a time, waiting for each to finish:
```
bash workflow/zea/run_all.sh
bash workflow/human/run_all.sh
bash workflow/monitoring/run_all.sh {prepare,chr7,panels,query,sweep,f0}
bash workflow/downsample/run_all.sh {prepare,bams,moth,models,query,pop,pop_call}
```
Notes:
- With no argument, `workflow/zea/run_all.sh` and `workflow/human/run_all.sh` run all their stages in order. For *H. zea* these are `reference align speedseq prepare sample ploidy model call call_y200`.
- `call_y200` repeats the Delly2 merges with Delly2's default minimum site quality (`-y 200`) for one comparison in the supplement. Every other Delly2 merge uses `-y 0`.
- The monitoring stage `f0` reruns every CNVcaller monitoring run with `-f 0` instead of `-f 0.1`, into a second results tree (`CNVB_RESULTS_F0`, default `results_f0/`), for Fig. S15.
- Jobs run on the local machine by default. To send each job to a scheduler, set `SUBMIT`, for example:
  ```
  SUBMIT="sbatch -c 8 --mem 32G -t 24:00:00 --wrap" bash workflow/zea/run_all.sh
  ```

**4. Make the tables and figures.** In the `cnvbench-analysis` environment, from the repository root:
```
bash analysis/make_all.sh
```
This runs every script in `analysis/` in the right order (some read the tables of others). Each script can also be run by itself, e.g. `python analysis/ddpcr.py`, and its docstring says what it makes. `trio_benchmark.py` caches its results in `results/tables/trio_benchmark.json`; after changing any calls, run it with `--recompute`.

## Where each figure and table comes from

Numbers follow the manuscript version of 22 November 2026. Supplementary items are in Additional file 1. Figures are written as PDF and PNG to `results/figures/`, and tables as TSV to `results/tables/`.

### Main figures

| Figure | Script | Output file |
|---|---|---|
| 1 Study design | `overview_figures.py` | `study_design` |
| 2 MIR through the filters, and by call length | `trio_benchmark.py` | `fig_filters_size_dup` |
| 3 Copy number at *Tryp77* against ddPCR | `ddpcr.py` | `fig_ddpcr_panels` |
| 4 Carrier-dilution panels | `monitoring_panels.py` | `timeseries_titration_nested_pooled_y0` |
| 5 Chromosome 7 monitoring scan | `monitoring_scan.py` | `monitoring_fp_headline` |

### Supplementary figures

| Figure | Script | Output file |
|---|---|---|
| S1 Match rule and filters (deletions) | `trio_benchmark.py` | `fig_metric_filter_del` |
| S2 Size bins and caller size ranges (deletions) | `trio_benchmark.py` | `fig_size_bins_del` |
| S3 MIR by call length (deletions) | `trio_benchmark.py` | `fig_mir_by_size_del` |
| S4 MIR with different background populations (deletions) | `trio_benchmark.py` | `fig_popdesign_del` |
| S5 Robustness of MIR to filters and match rules (deletions) | `trio_benchmark.py` | `fig_robustness_del` |
| S6 Caller consistency in *H. zea* and human trios (deletions) | `trio_benchmark.py` | `fig_species_del` |
| S7 *H. zea* parental calls longer than 1 Mb | `trio_benchmark.py` | `fig_oversized_calls` |
| S8–S11 Read depth, GC content and mapping quality around *Tryp77*, Trio1–Trio4 | `coverage_plots.py` | `Trio1_plot` … `Trio4_plot` |
| S12 CNVcaller frequency sweep | `monitoring_sweep.py` | `fig_freq_detect_published` |
| S13 Trio moths against the carrier-dilution panels | `monitoring_query.py` | `gatk_carrier_titration_nested` |
| S14 Trio moths against panels already held | `monitoring_query.py` | `fixed_query_stacked_y0` |
| S15 CNVcaller with and without its frequency threshold | `cnvcaller_f0.py` | `cnvcaller_f0_flipping` |
| S16 Size bins and caller size ranges (duplications) | `trio_benchmark.py` | `fig_size_bins_dup` |
| S17 Monitoring panels with other year mixes | `monitoring_panels.py` | `timeseries_monitoring_permutations_y0` |
| S18 Monitoring scan with three ways of running Delly2 | `monitoring_scan.py` | `monitoring_fp_delly` |
| S19 Calling at field depth (chromosome 7) | `downsample_chr7.py` | `downsample_chr7` |
| S20 MIR through the filters, and by call length (deletions) | `trio_benchmark.py` | `fig_filters_size_del` |
| S21 Match rule and filters (duplications) | `trio_benchmark.py` | `fig_metric_filter_dup` |
| S22 MIR with different background populations (duplications) | `trio_benchmark.py` | `fig_popdesign_dup` |
| S23 Robustness of MIR to filters and match rules (duplications) | `trio_benchmark.py` | `fig_robustness_dup` |
| S24 Caller consistency in *H. zea* and human trios (duplications) | `trio_benchmark.py` | `fig_species_dup` |
| S25 Read depth and duplication calls around *Tryp77* | `ddpcr.py` | `locus_depth_calls` |
| S26 Monitoring panel growing year by year | `monitoring_panels.py` | `timeseries_monitoring_sequential_y0` |

### Supplementary tables

| Table | Script | Output file |
|---|---|---|
| S1 Callers | (description, no script) | |
| S2 Deletion MIR under other filters | `trio_benchmark.py` | `mir_robustness` (deletion rows) |
| S3 MIR of each trio | `trio_benchmark.py` | `mir_headline` |
| S4 Callers compared within trios | `trio_benchmark.py` | `mir_paired` |
| S5 Human calls against HGSVC2 | `human_truthset.py` | `truthset_human` (HGSVC2 rows) |
| S6 Human deletions against HGSVC3 and Platinum Pedigree | `human_truthset.py` | `truthset_human` (other rows) |
| S7 CNVkit protocols | `trio_benchmark.py` | `mir_cnvkit_protocol` |
| S8 MIR without parental calls over 1 Mb | `trio_benchmark.py` | `mir_oversized` |
| S9 Chromosome map between assemblies | `chromosome_map.py` | `chromosome_map` |
| S10 ddPCR copy number | (measured, no script) | [`data/ddpcr_tryp77.tsv`](data/ddpcr_tryp77.tsv) |
| S11 Published *H. zea* sequencing data | `datasets.py` | `datasets` |
| S12 Published Louisiana moths | `monitoring_panels.py` | `monitoring_moths` |
| S13 Frequency-sweep populations | `monitoring_sweep.py` | `monitoring_sweep_populations` |
| S14 Carrier-dilution panels | `monitoring_panels.py` | `monitoring_nested_panels` |
| S15 Trio moths against panels already held | `monitoring_query.py` | `monitoring_fixed_query` |
| S16 Sequencing and alignment quality | `data_quality.py` | `data_quality` |
| S17 *Tryp77* copy number of each trio moth | `ddpcr.py` | `tryp77_cn` |
| S18 Duplication MIR under other filters | `trio_benchmark.py` | `mir_robustness` (duplication rows) |
| S19 Run time of each caller's main steps | `runtime.py` | `runtime_fixed` |
| S20 Run time and memory | `runtime.py` | `runtime_fixed`, `runtime` |
| S21 Chromosome 7 depth series | `downsample_chr7.py` | `downsample_chr7_summary` |
| S22 Drift-adjusted monitoring scan | `monitoring_scan.py` | `monitoring_scan_drift` |
| S23 Unrelated-parent control | `unrelated_parents.py` | `unrelated_parents` |
| S24 *Tryp77* copy number by background | `ddpcr.py` | `tryp77_ddpcr_by_background` |
| S25 CNVkit duplications against read depth | `cnvkit_depth.py` | `cnvkit_dups_depth` (also as `.tex`) |

### Numbers quoted in the text

| Script | Output file | What it gives |
|---|---|---|
| `trio_benchmark.py` | `mir_headline`, `filter_stages` | pooled MIR, and calls kept at each filter stage |
| `ddpcr.py` | `tryp77_ddpcr_correlation` | correlation of each caller with ddPCR |
| `tryp77_region.py` | `tryp77_region` | duplication calls around *Tryp77*, and the parental calls over 1 Mb (Section S3.3) |
| `blacklist_overlap.py` | `blacklist_overlap` | MIR at other blacklist overlaps (Methods) |
| `delly2_extras.py` | `delly2_extras` | per-sample Delly2, and `-y 0` against `-y 200` |
| `reference_stats.py` | `reference_stats` | share of each reference repeat-masked, gene lengths (Methods) |
| `monitoring_panels.py` | `monitoring_panel_calls`, `monitoring_nested_summary` | carriers called and false calls per panel |
| `monitoring_scan.py` | `monitoring_scan` | rank of *Tryp77* in the monitoring scan |
| `cnvcaller_f0.py` | `cnvcaller_f0_calls`, `cnvcaller_f0_summary` | CNVcaller at `-f 0` against `-f 0.1` (Section S4.4) |

Some scripts also write figures that are not in the paper: duplication versions of deletion-only supplementary figures, and `caller_scorecard` (from `overview_figures.py`).

## Run times

The run times of the trio benchmark come from re-timing every step at one configuration ([`data/runtime/fixed_cores.tsv`](data/runtime/fixed_cores.tsv)). Each step ran three times on 8 cores of one CPU type, with `workflow/runtime/timed.sh`. Each run called the Trio1 child (*H. zea*) or NA19240 (human) against crosstrio-9.

The run times of the published-moth preparation and the monitoring analyses are the records of the runs reported in the paper ([`data/runtime/steps.tsv`](data/runtime/steps.tsv)). They come from GNU `time` logs, Slurm accounting and worker timestamps.

## Tools

| Tool | Version |
|---|---|
| Delly2 | 1.7.2 |
| smoove | 0.2.8 (container `brentp/smoove:v0.2.8`, with LUMPY and SVTyper inside) |
| SpeedSeq (BWA-MEM alignment) | 0.1.2 (commit 4e60002) |
| GATK | 4.6.2.0, run with `--max-copy-number 10` |
| CNVcaller | commit d851794, with `workflow/callers/cnvcaller.patch` |
| CNVkit | 0.9.13 |
| Bowtie2 | 2.5.4 |
| Trimmomatic | 0.38 (trio moths; also the 2002 moths' read-depth BAMs), 0.39 (other published-moth trimming) |
| Python (analysis) | 3.13, with NumPy 2.3.0, SciPy 1.17.0, matplotlib 3.10.3, pysam 0.23.3 |

## Citation and licence

If you use this code or the ddPCR data, please cite the paper (citation to follow on publication).

Licence: MIT (see [`LICENSE`](LICENSE)). A Zenodo archive with a DOI will be added with the first release.
