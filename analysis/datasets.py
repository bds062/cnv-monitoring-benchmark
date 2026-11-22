"""Additional file 1, Table S11 (published H. zea genomic sequencing data): individuals, read length and sequenced depth per
BioProject, from the SRA run info in data/sra_runinfo.tsv. Which runs make up each row is set in ROWS below:
genomic DNA of H. zea only (no RNA-seq, RAD-seq, amplicon or other library strategies) and, where a BioProject
holds unpublished submissions, only the published ones. The trio children (this study) are not yet in the SRA;
their row comes from the read-pair counts in data/trio_children_reads.tsv. Table datasets.
    python analysis/datasets.py            # summarise data/sra_runinfo.tsv
    python analysis/datasets.py --fetch    # first refresh data/sra_runinfo.tsv from NCBI E-utilities"""
import collections
import csv
import io
import sys
import time
import urllib.request

from cnvbench.paths import REPO, TABLES

RUNINFO = REPO / "data" / "sra_runinfo.tsv"
CHILDREN = REPO / "data" / "trio_children_reads.tsv"
GENOME = 375e6            # bp; sequenced depth = sequenced bases per individual (or pool) / GENOME
COLS = ["Run", "BioProject", "BioSample", "SampleName", "LibraryName", "LibraryStrategy", "LibrarySource",
        "LibraryLayout", "Platform", "Model", "ScientificName", "ReleaseDate", "spots", "bases", "avgLength"]

# (section, label, BioProjects, run filter, field naming the individual or pool)
ROWS = [
    ("individuals", "PRJNA599999 (Fritz et al. 2020)", ["PRJNA599999"], lambda r: r["LibraryStrategy"] == "WGS", "BioSample"),
    # one BioSample per collection; each run (library) is one moth
    ("individuals", "PRJNA751583 (Taylor et al. 2021)", ["PRJNA751583"], lambda r: r["LibraryStrategy"] == "WGS", "LibraryName"),
    ("individuals", "PRJNA825115 (North et al. 2024)", ["PRJNA825115"], lambda r: r["LibraryStrategy"] == "WGS", "BioSample"),
    ("individuals", "PRJNA1055981 (Pezzini et al. 2024)", ["PRJNA1055981"], lambda r: r["LibraryStrategy"] == "WGS", "BioSample"),
    ("individuals", "PRJNA1126687 (Taylor et al. 2024)", ["PRJNA1126687"], lambda r: r["LibraryStrategy"] == "WGS", "BioSample"),
    # the five published genomes (released June 2024); later H. zea submissions are unpublished
    ("individuals", "PRJNA973566 (Trujillo et al. 2024)", ["PRJNA973566"],
     lambda r: r["LibraryStrategy"] == "WGS" and r["ReleaseDate"].startswith("2024"), "BioSample"),
    ("pools", "PRJNA767434 (Benowitz et al. 2022)", ["PRJNA767434"],
     lambda r: r["LibraryStrategy"] == "WGS" and r["Platform"] == "ILLUMINA", "BioSample"),
    ("pools", "PRJNA1111145 (Legan et al. 2024)", ["PRJNA1111145"], lambda r: r["LibraryStrategy"] == "WGS", "BioSample"),
    ("assemblies", "PRJNA378438 (Pearce et al. 2017)", ["PRJNA378438"], lambda r: r["LibraryStrategy"] == "WGS", "BioSample"),
    ("assemblies", "PRJNA767434 (Benowitz et al. 2022), ASM2234304v1", ["PRJNA767434"],
     lambda r: r["LibraryStrategy"] == "WGS" and r["Platform"] == "PACBIO_SMRT", "BioSample"),
    ("assemblies", "PRJNA804956, PRJNA788876 (Stahlke et al. 2023), ilHelZeax1.1", ["PRJNA804956", "PRJNA788876"],
     lambda r: True, "BioSample"),
]
PLATFORM = {"PACBIO_SMRT": "PacBio", "LS454": "454"}


def fetch():
    """Refresh data/sra_runinfo.tsv from NCBI E-utilities (esearch + efetch rettype=runinfo)."""
    eu = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    rows = []
    for bp in sorted({b for row in ROWS for b in row[2]}):
        es = urllib.request.urlopen(f"{eu}esearch.fcgi?db=sra&term={bp}%5BBioProject%5D&retmax=5000&usehistory=y").read().decode()
        we, qk = es.split("<WebEnv>")[1].split("<")[0], es.split("<QueryKey>")[1].split("<")[0]
        txt = urllib.request.urlopen(f"{eu}efetch.fcgi?db=sra&query_key={qk}&WebEnv={we}&rettype=runinfo&retmode=text").read().decode()
        rows += [{c: r[c] for c in COLS} for r in csv.DictReader(io.StringIO(txt)) if r.get("Run", "")[:3] in ("SRR", "ERR", "DRR")]
        time.sleep(0.5)
    rows.sort(key=lambda r: (r["BioProject"], r["Run"]))
    with open(RUNINFO, "w") as fh:
        fh.write(f"# NCBI SRA run info (E-utilities efetch rettype=runinfo) for the BioProjects of Table S11, fetched "
                 f"{time.strftime('%Y-%m-%d')}.\n# Refresh: python analysis/datasets.py --fetch\n")
        w = csv.DictWriter(fh, COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("wrote", RUNINFO, len(rows), "runs")


def rng(vals, fmt="{:.0f}"):
    lo, hi = fmt.format(min(vals)), fmt.format(max(vals))
    return lo if lo == hi else f"{lo}-{hi}"


def read_length(runs):
    out = []
    short = [r for r in runs if r["Platform"] not in PLATFORM]
    for paired in (True, False):
        rs = [r for r in short if (r["LibraryLayout"] == "PAIRED") == paired]
        if rs:
            mate = [float(r["avgLength"]) / (2 if paired else 1) for r in rs]
            out.append(("2x" if paired else "") + rng(mate))
    for plat, name in PLATFORM.items():
        rs = [r for r in runs if r["Platform"] == plat]
        if rs:                                              # mean read length over the platform's bases
            mean = sum(int(r["bases"]) for r in rs) / sum(int(r["spots"]) for r in rs)
            out.append(f"{name} {mean / 1e3:.1f} kb")
    return "; ".join(out)


def main():
    if "--fetch" in sys.argv:
        fetch()
    runs = [r for r in csv.DictReader((l for l in open(RUNINFO) if not l.startswith("#")), delimiter="\t")
            if r["ScientificName"] == "Helicoverpa zea"]
    out = []
    for section, label, bps, keep, unit in ROWS:
        rs = [r for r in runs if r["BioProject"] in bps and keep(r)]
        bases = collections.Counter()
        for r in rs:
            bases[r[unit]] += int(r["bases"])
        depth = [b / GENOME for b in bases.values()]
        row = dict(section=section, bioproject=label, n=len(bases), runs=len(rs), read_length_bp=read_length(rs),
                   depth=rng(depth) + "x")
        out.append(row)
        print(f"{section:11s} {label:58s} n={row['n']:<4d} runs={row['runs']:<4d} {row['read_length_bp']:24s} {row['depth']}")
    # trio children (this study): read pairs from data/trio_children_reads.tsv until the SRA release
    kids = list(csv.DictReader((l for l in open(CHILDREN) if not l.startswith("#")), delimiter="\t"))
    depth = [2 * int(k["read_pairs"]) * int(k["read_length"]) / GENOME for k in kids]
    row = dict(section="individuals", bioproject="This study (trio children)", n=len(kids), runs=len(kids),
               read_length_bp="2x" + rng([int(k["read_length"]) for k in kids]), depth=rng(depth) + "x")
    out.insert(6, row)
    print(f"{row['section']:11s} {row['bioproject']:58s} n={row['n']:<4d} runs={row['runs']:<4d} "
          f"{row['read_length_bp']:24s} {row['depth']}")
    TABLES.mkdir(parents=True, exist_ok=True)
    with open(TABLES / "datasets.tsv", "w") as fh:
        w = csv.DictWriter(fh, list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    print("wrote", TABLES / "datasets.tsv")


if __name__ == "__main__":
    main()
