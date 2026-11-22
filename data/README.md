# Data

This folder holds two things:
- `download.sh`, which fetches the large public inputs: references, reads, truth sets and tools;
- small tables that are kept in the repository: the ddPCR copy numbers, the sample sheets and a few reference tables.

## ddPCR copy number at *Tryp77* (`ddpcr_tryp77.tsv`)

[`ddpcr_tryp77.tsv`](ddpcr_tryp77.tsv) gives the copy number of *Tryp77* in each of the 12 *H. zea* trio moths, measured by droplet digital PCR (ddPCR). It is Table S10 of the paper (Additional file 1), and every analysis that compares callers with ddPCR reads this file.

**Assay:**
- Instrument: Bio-Rad QX droplet digital PCR.
- Probes: a FAM probe on an amplicon in *Tryp77*, and a HEX probe on the single-copy reference gene *HelRegB*.
- The amplicon lies at NC_061458.1:8,252,593–8,252,690 of the reference GCF_022581195.2.
- Copy number = 2 × (*Tryp77* / *HelRegB*), the `CNV` column of the instrument export.

| Column | Meaning |
|---|---|
| `trio` | Trio1–Trio4 |
| `role` | `resistant` parent, `susceptible` parent or `child` |
| `sample` | Name of the moth's sequencing library (as in `samples/trios.tsv`) |
| `ddpcr_label` | Name of the sample on the ddPCR plate |
| `run_date` | Date of the ddPCR run (YYYY-MM-DD) |
| `well` | Plate well |
| `ratio` | *Tryp77* / *HelRegB* ratio from the instrument export |
| `cn` | *Tryp77* copy number, 2 × `ratio` (rounded to two decimals) |

**Notes:**
- The ddPCR plate labels were matched to the sequenced moths by hand. Each child is the "_2" sibling of its cross; for example, ddPCR sample BV20_AA27_P6BZM_C3_2 is sequenced child BV20_AA27_P6BZM_C32.
- The Trio2 child has no ddPCR well (its columns read `NA`). Its copy number of 2 comes from its pedigree and from read-depth genotyping of its F2 offspring (see the paper's Methods).
- The analysis uses 2 × `ratio` at full precision, and `cn` where there is no ratio.

## Inputs fetched by `download.sh`

All inputs are public. `bash data/download.sh [TARGET ...]` fetches them into `$CNVB_DATA` (default `data/`, see `config.sh`). With no target it fetches everything. The read downloads need the environment `envs/preprocess.yml`.

| Target | What | Source | Saved to |
|---|---|---|---|
| `zea_ref` | *H. zea* reference ilHelZeax1.1 (soft-masked by NCBI), and the earlier GA-R assembly (for the trio coverage, GC and mapping-quality comparison) | NCBI RefSeq GCF_022581195.2, GenBank GCA_022343045.1 | `ref/zea/`, `ref/zea_gca/` |
| `trios` | Reads of the 12 trio moths | Parents: PRJNA1126687 (Taylor et al. 2024). Children: this study (released in the SRA with the paper) | `reads/zea/{sample}_{1,2}.fastq.gz` |
| `published` | 35 Louisiana moths (2002, 2012, 2017) and 6 Benzon cross founders | PRJNA751583 (Taylor et al. 2021), PRJNA1126687 | `reads/zea/{SRR}_{1,2}.fastq.gz` |
| `human_ref` | GRCh38 analysis set used by the 1000 Genomes alignments | `GRCh38_full_analysis_set_plus_decoy_hla.fa` | `ref/human/` |
| `crams` | 12 high-coverage 1000 Genomes CRAMs (md5-checked) | ENA; Byrska-Bishop et al. 2022 | `reads/human/{sample}.final.cram` |
| `truth` | Long-read SV call sets for the human children | HGSVC2 freeze 4, HGSVC3 v1.0, Platinum Pedigree v1.2 (NA12878) | `truth/{hgsvc2,hgsvc3,platinum}/` |
| `repeats` | UCSC RepeatMasker track for hg38 | `rmsk.txt.gz` | `repeats/hg38_rmsk.txt.gz` |
| `tools` | Trimmomatic 0.38 and 0.39; GATK 4.6.2.0; CNVcaller d851794 (patched with `workflow/callers/cnvcaller.patch`); SpeedSeq 0.1.2 (commit 4e60002); the smoove 0.2.8 image | GitHub, Docker Hub | `tools/` |

After the `tools` target, two manual steps are needed:
1. Install GATK's Python package into the `gatk` environment:
   `pip install $CNVB_DATA/tools/gatk-4.6.2.0/gatkPythonPackageArchive.zip`
2. Put `$CNVB_DATA/tools/speedseq/bin` on `PATH` in the `speedseq` environment.

Tool locations can be overridden with:
- `TRIMMOMATIC_DIR` (the folder holding `Trimmomatic-0.38/` and `Trimmomatic-0.39/`);
- `GATK`;
- `CNVCALLER_DIR`;
- `SMOOVE_SIF`.

The repeat BED files used by the MIR repeat filter (`repeats/zea_softmask.bed`, `repeats/human_rmsk.bed`) are made by `workflow/{zea,human}/prepare_reference.sh`.

## Sample sheets (`samples/`)

- **`trios.tsv`**: species, trio, role and sample name of the 24 trio members (12 moths, 12 people).
- **`trio_runs.tsv`**: SRA run of each trio moth.
  - The eight parents are cross founders published by Taylor et al. (2024) under their lab names.
  - The four children are released with this study. Until their accessions exist, their rows read `pending` and `download.sh` skips them.
- **`published_moths.tsv`**: the 41 published moths. Columns: sample (= SRA run), experiment, BioProject, group (LA2002, LA2012, LA2017, Benzon), and the trimming of the moth's read-depth alignment, which `workflow/zea/align.sh` sets by group:
  - `trio` for the 2002 moths: the trio moths' settings.
  - `taylor2021` for the 2012 and 2017 moths: Trimmomatic 0.39 with Taylor et al.'s (2021) 34-sequence TruSeq file, `workflow/zea/adapt_seq.fa`. Trimmomatic uses only 2 of its sequences and no palindrome pair.
  - `taylor2024` for the Benzon founders: Trimmomatic 0.39, 4 bp windows at quality 30, reads of at least 36 bp. `align.sh` fills in the values that paper does not give.
  - The paper's founder runs used the reads as Taylor et al. trimmed them. Those trimmed reads are not public, so this workflow re-trims the public reads instead.
  - The 2002 moths and the founders also get a second Bowtie2 alignment with the `taylor2021` trimming (`align.sh SAMPLE taylor2021`, written to `$CNVB_RESULTS/t21`). The monitoring analyses use it for the 2002 moths, and the BWA-MEM realignments for Delly2 and smoove are made from it.
- **`backgrounds.tsv`**: the members of each background population. Columns: species, background, query trio (`all` if every sample shares it), member.
  - `crosstrio9`: the nine moths or people of the other three trios.
  - `crosstrio3`: one parent from each other trio.
    - *H. zea*: the susceptible parents.
    - Human: the parent whose sex chromosomes match the *H. zea* parent in that trio slot (YRI mother, PJL father, CEU mother, PUR father).
  - `crosstrio3_auto` (*H. zea*, CNVkit only): the crosstrio3 moths under the earlier CNVkit protocol (automatic 4,589 bp bins).
  - `pooled12`: all 12 trio moths.
  - `susceptible13`: the 13 Louisiana 2002 moths.
  - `neutral18` (balanced-18 in the paper), 18 moths:
    - six Louisiana 2002 moths (SRX11654257–262);
    - six Louisiana 2017 moths (SRX11764981–986);
    - six Benzon founders, one from each Taylor et al. (2024) cross not used as a trio (libraries IL100235829, 830, 832, 835, 838 and 847).
  - `neutral9` (balanced-9): the first three of each balanced-18 group.
- **`human_crams.tsv`**: CRAM URL and md5 checksum of each human sample.
- **`timeseries_panels.tsv`**: the members of every monitoring panel: the year panels, the nested carrier-dilution panels (`nest_*`) and the frequency-sweep populations (`sweep_*`).

### Benzon founders: library names vs. SRA names

The Benzon founders were sequenced as libraries IL100235828–847 but deposited under cross names. Matching by shared reads gives:

| Library | SRA sample | Run |
|---|---|---|
| IL100235829 | BV20_AA27_P8_BZM2 | SRR29490991 |
| IL100235830 | BV20_BP68_O8_BZM1 | SRR29490980 |
| IL100235832 | BV20_AA27_O8_BZF1 | SRR29490978 |
| IL100235835 | BV20_G7_04_BZF2 | SRR29490975 |
| IL100235838 | BV20_AA27_O1_BZF1 | SRR29490990 |
| IL100235847 | BV_DB104_O6_BZF2 | SRR29490981 |

In each founder pair the SRA `breed` labels are swapped: the moth carrying the *Tryp77* amplification is labelled susceptible.

## Other tables

| File | What | Used by |
|---|---|---|
| `sra_runinfo.tsv` | NCBI SRA run information for every published *H. zea* genomic BioProject (refresh with `python analysis/datasets.py --fetch`) | `datasets.py` (Table S11) |
| `trio_children_reads.tsv` | Sequenced read pairs of the four trio children, until they are in the SRA | `datasets.py` (Table S11) |
| `ogs2_gene_lengths.tsv` | Lengths of the *H. zea* OGS2 gene models, lifted onto GCF_022581195.2; they set the MIR size bins | `reference_stats.py`, `trio_benchmark.py` |
| `wfmash/alignments.tsv.gz`, `wfmash/chromosome_names.tsv` | wfmash alignments of the earlier GA-R assembly to the reference, and the chromosome name of each accession | `chromosome_map.py` (Table S9) |
| `runtime/fixed_cores.tsv` | Re-timing of every trio-benchmark step: 8 cores, three repeats | `runtime.py` (Tables S19, S20) |
| `runtime/steps.tsv` | Run times recorded in the paper's runs of the published-moth preparation and the monitoring analyses | `runtime.py` (Table S20) |

`sra_runinfo.tsv`, `trio_children_reads.tsv`, `ogs2_gene_lengths.tsv` and the `runtime/` tables start with comment lines (`#`) that describe their columns and source. The sample sheets and the ddPCR table have no comment lines, so GitHub shows them as tables; they are described above.
