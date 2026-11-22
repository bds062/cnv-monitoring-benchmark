"""Trio sample sheet and one loader per caller output format.

Every loader returns Call tuples on normalised chromosome names ('1'..'31' in H. zea, with the
Z chromosome as '31'; '1'..'22', X = '23', Y = '24' in human): deletions and duplications of at
least MIN_LEN bp with a non-reference genotype. VCFs may be plain or gzipped."""
import csv
import gzip
from collections import namedtuple

from .paths import SAMPLES, call_dir

Call = namedtuple("Call", "chrom start end svtype")
Trio = namedtuple("Trio", "name parents child")       # parents: (resistant, susceptible) or (father, mother)
CALLERS = ["Delly2", "smoove", "GATK-gCNV", "CNVcaller", "CNVkit"]     # breakpoint | read depth
DIR = {"Delly2": "delly2", "smoove": "smoove", "GATK-gCNV": "gatk", "CNVcaller": "cnvcaller", "CNVkit": "cnvkit"}
HEADLINE = {"Delly2": "crosstrio9", "smoove": "crosstrio9", "GATK-gCNV": "crosstrio9", "CNVcaller": "crosstrio9",
            "CNVkit": "crosstrio3"}
VCF = {"delly2": "{}.vcf", "smoove": "{}.vcf", "gatk": "{}.segments.vcf", "cnvcaller": "{}.vcf",
       "cnvkit": "{}.call.vcf", "delly2_cnv": "{}.vcf"}
MIN_LEN = 100
SPECIES = {"zea": "H. zea", "human": "Human"}


def trios(species):
    rows = [r for r in csv.DictReader(open(SAMPLES / "trios.tsv"), delimiter="\t") if r["species"] == species]
    out = []
    for name in dict.fromkeys(r["trio"] for r in rows):
        by = {r["role"]: r["sample"] for r in rows if r["trio"] == name}
        par = ("resistant", "susceptible") if species == "zea" else ("father", "mother")
        out.append(Trio(name, (by[par[0]], by[par[1]]), by["child"]))
    return out


def samples(species):
    return [s for t in trios(species) for s in (*t.parents, t.child)]


def roles(species):
    return {r["sample"]: r["role"] for r in csv.DictReader(open(SAMPLES / "trios.tsv"), delimiter="\t")
            if r["species"] == species}


def norm_chrom(species, name):
    """Normalised chromosome name, or None for scaffolds, organelles and alt contigs."""
    name = name[3:] if name.startswith("chr") else name
    if species == "zea":
        if name in ("X", "Z"):                     # CNVkit is given the Z chromosome as X
            return "31"
        if name.startswith("NC_0614"):             # RefSeq NC_061452.1 .. NC_061482.1 = chromosomes 1-31
            n = int(name[3:9]) - 61451
        elif name.isdigit():
            n = int(name)
        else:
            return None
        return str(n) if 1 <= n <= 31 else None
    if name in ("X", "Y"):
        return "23" if name == "X" else "24"
    return name if name.isdigit() and 1 <= int(name) <= 22 else None


def open_text(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def find(folder, name):
    """folder/name.gz if present, else folder/name."""
    gz = folder / (name + ".gz")
    return gz if gz.exists() else folder / name


def vcf_path(species, caller, background, sample):
    d = DIR.get(caller, caller)
    return find(call_dir(species, d, background, sample), VCF[d].format(sample))


def vcf_records(path):
    """(header columns, record fields) for every record."""
    cols = None
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                cols = line[1:].rstrip("\n").split("\t")
                continue
            yield cols, line.rstrip("\n").split("\t")


def info(field):
    return dict(kv.split("=", 1) for kv in field.split(";") if "=" in kv)


def fmt(rec, col=9):
    return dict(zip(rec[8].split(":"), rec[col].split(":")))


def nonref(gt):
    return gt not in ("./.", "0/0", ".", ".|.", "0|0")


def sample_column(cols, sample):
    """CNVcaller names its columns after the input files, so match the sample name or a prefix of it."""
    return next(i for i, c in enumerate(cols) if i >= 9 and (c == sample or c.startswith(sample + "_")))


# ─── one function per output format: (cols, rec, INFO, sample) -> call types of this record ─────
def _breakpoint(cols, r, inf, sample):            # delly call and smoove (LUMPY + SVTyper); one sample per file
    return (inf.get("SVTYPE"),) if nonref(fmt(r)["GT"]) else ()


def _delly2(cols, r, inf, sample):
    return _breakpoint(cols, r, inf, sample) if r[6] == "PASS" else ()


def _gatk(cols, r, inf, sample):
    # GATK-gCNV writes GT ./. for duplication segments, so gains are read from CN
    f = fmt(r)
    if r[4] == "<DEL>":
        return ("DEL",) if nonref(f["GT"]) else ()
    if r[4] == "<DUP>":
        return ("DUP",) if f.get("CN", ".") not in (".", "") and int(f["CN"]) > 2 else ()
    return ()


def _cnvcaller(cols, r, inf, sample):
    # One record per CNV region; each individual's call is typed by its own alleles (CN0 = deletion,
    # CN2/CNH = gain), so one individual can carry a deletion and a duplication at the same region.
    if inf.get("SVTYPE") not in ("DEL", "DUP"):
        return ()
    gt = fmt(r, sample_column(cols, sample))["GT"]
    if not nonref(gt):
        return ()
    alts = r[4].split(",")
    idx = [int(a) for a in gt.replace("|", "/").split("/") if a not in (".", "0")]
    return sorted({"DEL" if alts[i - 1] == "CN0" else "DUP" for i in idx})


PARSE = {"delly2": _delly2, "smoove": _breakpoint, "gatk": _gatk, "cnvcaller": _cnvcaller,
         "cnvkit": _breakpoint}


def diploid_z(species, caller, background):
    """GATK-gCNV and the earlier CNVkit protocol report the single Z chromosome of ZW females as a deletion
    of most of the chromosome, so their H. zea Z calls are left out (Methods)."""
    return species == "zea" and (caller == "GATK-gCNV" or (caller == "CNVkit" and background != "crosstrio3"))


def load(species, caller, background, sample, min_len=MIN_LEN):
    parse = PARSE[DIR.get(caller, caller)]
    skip_z = diploid_z(species, caller, background)
    out = []
    for cols, r in vcf_records(vcf_path(species, caller, background, sample)):
        chrom = norm_chrom(species, r[0])
        if chrom is None or (skip_z and chrom == "31"):
            continue
        inf = info(r[7])
        start, end = int(r[1]), int(inf.get("END", r[1]))
        if end - start + 1 < min_len:
            continue
        out += [Call(chrom, start, end, t) for t in parse(cols, r, inf, sample) if t in ("DEL", "DUP")]
    return out


def load_set(species, designs):
    """{label: (caller, background)} -> {label: {sample: [Call]}} for all trio samples."""
    return {lab: {s: load(species, c, bg, s) for s in samples(species)} for lab, (c, bg) in designs.items()}


def headline(species):
    return load_set(species, {c: (c, HEADLINE[c]) for c in CALLERS})
