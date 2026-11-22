"""The Tryp77 locus on H. zea chromosome 7: ddPCR copy number, each caller's copy number over the ddPCR
amplicon, Delly2 read-depth windows and read depth from BAMs."""
import csv
import math
import tempfile

import pysam

from .calls import find, fmt, info, nonref, norm_chrom, open_text, sample_column, vcf_path, vcf_records
from .paths import REPO, call_dir

CHROM = "7"                                     # NC_061458.1
AMPLICON = (8_252_593, 8_252_690)               # ddPCR primer span
GENES = [("Tryp079", 8_229_974, 8_233_384), ("Tryp078", 8_241_790, 8_246_612), ("Tryp077", 8_250_592, 8_254_426)]
TRYP = (GENES[0][1], GENES[-1][2])               # Tryp079-Tryp077


def ddpcr():
    """{sample: Tryp77 copy number} for the 12 trio moths (data/ddpcr_tryp77.tsv): 2 x the Tryp77/HelRegB ratio at
    full precision, or the table's copy number where there is no ddPCR well (the Trio2 child, ratio NA)."""
    rows = csv.DictReader(open(REPO / "data" / "ddpcr_tryp77.tsv"), delimiter="\t")
    return {r["sample"]: 2 * float(r["ratio"]) if r["ratio"] != "NA" else float(r["cn"]) for r in rows}


def most_extreme(events, a, b, bin_size=200):
    """Copy number over [a, b): in each 200 bp bin the overlapping event furthest from two copies (2 if none),
    averaged over the bins."""
    vals = []
    for i in range(math.ceil((b - a) / bin_size)):
        s, e = a + i * bin_size, min(a + (i + 1) * bin_size, b)
        best = 2.0
        for x, y, cn in events:
            if y > s and x < e and abs(cn - 2) > abs(best - 2):
                best = cn
        vals.append(best)
    return sum(vals) / len(vals)


def _overlapping(species, caller, background, sample, chrom, a, b):
    for cols, r in vcf_records(vcf_path(species, caller, background, sample)):
        if norm_chrom(species, r[0]) == chrom:
            end = int(info(r[7]).get("END", r[1]))
            if end > a and int(r[1]) < b:
                yield cols, r, end


def caller_cn(caller, background, sample, a=AMPLICON[0], b=AMPLICON[1], chrom=CHROM, species="zea"):
    """Copy number one caller reports for one moth over [a, b) of a chromosome (default: the ddPCR amplicon)."""
    if caller == "Delly2" and background == "none":     # read-depth copy number of `delly cnv`, one moth
        num = den = 0.0
        for c, s, e, cn in depth_windows(sample, species):
            if c == chrom and e > a and s < b:
                L = min(e, b) - max(s, a)
                num, den = num + cn * L, den + L
        return num / den if den else math.nan
    if caller == "CNVcaller":                   # continuous copy number of the CNV region over the amplicon
        for cols, r, _ in _overlapping(species, caller, background, sample, chrom, a, b):
            return float(fmt(r, sample_column(cols, sample))["CP"])
        return 2.0
    if caller == "Delly2":                      # population `delly cnv` workflow against a background
        return delly_pop_cn(background, sample, a, b, chrom, species)
    ev = []
    for cols, r, end in _overlapping(species, caller, background, sample, chrom, a, b):
        f = fmt(r)
        if caller == "smoove":                # no copy number: two-allele genotype as a proxy
            sv = info(r[7]).get("SVTYPE")
            if sv in ("DEL", "DUP") and nonref(f["GT"]):
                hom = f["GT"] in ("1/1", "1|1")
                ev.append((int(r[1]), end, {"DEL": (1, 0), "DUP": (3, 4)}[sv][hom]))
        elif f.get("CN", ".") != "." and (caller == "GATK-gCNV" or r[4] in ("<DEL>", "<DUP>")):
            ev.append((int(r[1]), end, float(f["CN"])))
    return most_extreme(ev, a, b)


def info_end(r):
    """INFO/END of a record: pysam's r.stop follows the BCF length field, which is 1 in these Delly2 BCFs."""
    info = str(r).split("\t")[7]
    return next(int(kv[4:]) for kv in info.split(";") if kv.startswith("END="))


def delly_pop_cn(background, sample, a=AMPLICON[0], b=AMPLICON[1], chrom=CHROM, species="zea"):
    """Integer CN of a moth at the shortest PASS site over [a, b) of Delly2's population read-depth workflow
    (workflow/callers/delly2_pop.sh: delly cnv -l, delly merge -y 0, delly cnv -u, delly classify -f germline),
    run on the moth with a background: delly2_pop/{background}/{sample}/{sample}.bcf; 2 if no site."""
    path = find(call_dir(species, "delly2_pop", background, sample), f"{sample}.bcf")
    best = None
    with pysam.VariantFile(str(path)) as vf:
        for r in vf:
            if norm_chrom(species, r.chrom) != chrom or "PASS" not in r.filter.keys():
                continue
            end = info_end(r)
            if end > a and r.pos < b and (best is None or end - r.pos < best[0]):
                best = (end - r.pos, r.samples[sample]["CN"])
    return 2.0 if best is None or best[1] is None else float(best[1])


_windows = {}


def depth_windows(sample, species="zea"):
    """[(chrom, start, end, CN)] read-depth windows of `delly cnv` (delly2_cnv/none/{sample}/{sample}.cov.gz)."""
    if sample not in _windows:
        rows = []
        with open_text(find(call_dir(species, "delly2_cnv", "none", sample), f"{sample}.cov")) as fh:
            col = [i for i, h in enumerate(fh.readline().rstrip("\n").split("\t")) if h.endswith("_CN")][-1]
            for line in fh:
                p = line.rstrip("\n").split("\t")
                c = norm_chrom(species, p[0])
                if c is not None:
                    rows.append((c, int(p[1]), int(p[2]), float(p[col])))
        _windows[sample] = rows
    return _windows[sample]


def bam_depth(bam, contig, windows):
    """Mean read depth in each (start, end) window (samtools bedcov)."""
    with tempfile.NamedTemporaryFile("w", suffix=".bed") as fh:
        fh.writelines(f"{contig}\t{a}\t{b}\n" for a, b in windows)
        fh.flush()
        out = pysam.bedcov(fh.name, str(bam))
    return [int(l.split("\t")[3]) / (int(l.split("\t")[2]) - int(l.split("\t")[1])) for l in out.strip().split("\n")]


def bam_contig(path, species, chrom):
    """Name of a normalised chromosome in a BAM header."""
    return next(c for c in pysam.AlignmentFile(str(path)).references if norm_chrom(species, c) == chrom)
