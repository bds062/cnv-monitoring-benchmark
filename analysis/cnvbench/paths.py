"""Repository paths. Data and results default to folders inside the repository;
set CNVB_DATA / CNVB_RESULTS to read or write them elsewhere (same variables as config.sh)."""
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = Path(os.environ.get("CNVB_DATA", REPO / "data"))
RESULTS = Path(os.environ.get("CNVB_RESULTS", REPO / "results"))
SAMPLES = REPO / "data" / "samples"          # sample sheets are tracked, so always in the repository
CALLS = RESULTS / "calls"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"


REF = {"zea": DATA / "ref" / "zea" / "GCF_022581195.2_ilHelZeax1.1_genomic.fna",
       "human": DATA / "ref" / "human" / "GRCh38_full_analysis_set_plus_decoy_hla.fa",
       "zea_gca": DATA / "ref" / "zea_gca" / "GCA_022343045.1_ASM2234304v1_genomic.fna"}   # earlier GA-R assembly
REPEATS = {"zea": DATA / "repeats" / "zea_softmask.bed", "human": DATA / "repeats" / "human_rmsk.bed"}
TRUTH = DATA / "truth"                       # hgsvc2/, hgsvc3/, platinum/ (human truth sets)


def call_dir(species, caller, background, sample):
    """results/calls/{species}/{caller}/{background}/{sample}/ (background 'none' for per-sample callers)."""
    return CALLS / species / caller / background / sample


def bam(species, sample, kind="bam"):
    """results/{kind}/{species}/{sample}.bam; kind = bam (Bowtie2 / 1000 Genomes), bam_bwa (SpeedSeq BWA-MEM)
    or bam_gca (Bowtie2 to the earlier GCA_022343045.1 assembly, H. zea only)."""
    return RESULTS / kind / species / f"{sample}.bam"


def trim_log(species, sample):
    """Trimmomatic log of one sample (H. zea trio moths)."""
    return RESULTS / "trim" / species / f"{sample}.trimmomatic.log"
