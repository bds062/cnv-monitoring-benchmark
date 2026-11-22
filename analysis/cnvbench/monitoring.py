"""Monitoring analyses: panels, the caller-free read-depth truth and each caller's copy number at the
Tryp77 ddPCR amplicon, read from what workflow/monitoring writes under results/monitoring/."""
import csv
import subprocess
from collections import defaultdict

from .calls import fmt, info, norm_chrom, open_text, vcf_records
from .locus import AMPLICON, CHROM, caller_cn
from .paths import RESULTS, SAMPLES, TABLES

MON = RESULTS / "monitoring"
TRYP77 = (8_219_000, 8_430_000)          # the ~200 kb read-depth CNV holding Tryp079-Tryp077 (delly cnv)
DUP, DEL = 2.5, 1.5                      # a call is a duplication at >= 2.5 copies, a deletion at <= 1.5
YEARS = ("2002", "2012", "2017")
CALLERS = ("GATK-gCNV", "CNVcaller", "Delly2 merge (-y 0)", "Delly2 merge (default)")
MINQUAL = {"Delly2 merge (-y 0)": 0, "Delly2 merge (default)": 200}


def _tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


PANELS = defaultdict(list)
for _r in _tsv(SAMPLES / "timeseries_panels.tsv"):
    PANELS[_r["panel"]].append(_r["member"])
YEAR = {s: y for y in YEARS for s in PANELS[y]}
TRIO_MOTHS = [(r["trio"], r["role"], r["sample"]) for r in _tsv(SAMPLES / "trios.tsv") if r["species"] == "zea"]


def read_depth():
    """{moth: (mean depth, read-depth copy number at the amplicon)} for the 35 published moths."""
    return {r["sample"]: (float(r["depth"]), float(r["depth_cn"])) for r in _tsv(MON / "read_depth.tsv")}


def end(rec):
    return int(info(rec[7])["END"])


def on_chr7(chrom):
    return norm_chrom("zea", chrom) == CHROM


def on_amplicon(chrom, start, stop):
    return on_chr7(chrom) and start <= AMPLICON[1] and stop >= AMPLICON[0]


def gatk_cn(path):
    """Integer CN of the GATK-gCNV segment over the amplicon."""
    return next(float(fmt(r)["CN"]) for _, r in vcf_records(path) if on_amplicon(r[0], int(r[1]), end(r)))


def cnvcaller_cp(path):
    """{member: CP} of the CNVcaller CNV region over the amplicon; 2 for everyone if there is none."""
    out = {}
    for cols, r in vcf_records(path):
        out = out or dict.fromkeys(cols[9:], 2.0)
        if on_amplicon(r[0], int(r[1]), end(r)):
            out = {n: float(fmt(r, 9 + i)["CP"]) for i, n in enumerate(cols[9:])}
    return out


def bcf(path, form, *opts):
    """bcftools query lines."""
    return subprocess.run(["bcftools", "query", *opts, "-f", form, str(path)], capture_output=True,
                          text=True, check=True).stdout.splitlines()


def delly_cn(filtered, genotypes=None):
    """{moth: integer CN} at the shortest PASS site over the amplicon after delly classify; 2 if none.
    genotypes: a moth genotyped at the panel's sites (fixed query) instead of the panel itself."""
    sites = [(int(b) - int(a), i) for c, i, a, b, f in (x.split("\t") for x in bcf(filtered, "%CHROM\t%ID\t%POS\t%INFO/END\t%FILTER\n"))
             if f == "PASS" and on_amplicon(c, int(a), int(b))]
    src = genotypes or filtered
    names = bcf(src, "", "-l")
    if not sites:
        return dict.fromkeys(names, 2.0)
    site = min(sites)[1]
    cns = next(x.split("\t")[1:] for x in bcf(src, "%ID[\t%CN]\n") if x.split("\t")[0] == site)
    return {n: 2.0 if c == "." else float(c) for n, c in zip(names, cns)}


def delly_alone(moth):
    """delly cnv copy number of a moth called alone on chromosome 7: length-weighted over the amplicon."""
    with open_text(MON / "chr7" / moth / f"{moth}.cov.gz") as fh:
        col = next(i for i, h in enumerate(fh.readline().rstrip("\n").split("\t")) if h.endswith("_CN"))
        w = [(min(int(p[2]), AMPLICON[1]) - max(int(p[1]), AMPLICON[0]), float(p[col]))
             for p in (line.split("\t") for line in fh)]
    w = [(L, c) for L, c in w if L > 0]
    return sum(L * c for L, c in w) / sum(L for L, _ in w)


def panel_cn(panel):
    """{caller: {member: copy number at the amplicon}} with every caller run on the panel itself."""
    d = MON / "panels" / panel
    out = {"GATK-gCNV": {s: gatk_cn(d / "gatk" / f"{s}.segments.vcf.gz") for s in PANELS[panel]},
           "CNVcaller": cnvcaller_cp(d / "cnvcaller" / f"{panel}.vcf.gz")}
    for c, y in MINQUAL.items():
        out[c] = delly_cn(d / f"delly_y{y}" / "filtered.bcf")
    return out


def query_cn(panel, moth, caller):
    """Copy number of a trio moth called against a panel (fixed query, carrier titration). crosstrio9: the
    trio-benchmark calls for GATK-gCNV and CNVcaller, and for Delly2 the merge workflow on the other trios."""
    if panel == "crosstrio9" and caller in ("GATK-gCNV", "CNVcaller"):
        return caller_cn(caller, "crosstrio9", moth)
    if caller == "GATK-gCNV":
        return gatk_cn(MON / "query" / panel / "gatk" / f"{moth}.segments.vcf.gz")
    if caller == "CNVcaller":
        return cnvcaller_cp(MON / "query" / panel / "cnvcaller" / f"{moth}.vcf.gz")[moth]
    if panel == "crosstrio9":
        panel = "crosstrio9_" + next(t for t, _, s in TRIO_MOTHS if s == moth)
    y = MINQUAL[caller]
    return delly_cn(MON / "panels" / panel / f"delly_y{y}" / "filtered.bcf",
                    MON / "query" / panel / f"delly_y{y}" / f"{moth}.bcf")[moth]


def outcome(carrier, cn):
    """tp/fn for a carrier, tn/fp (false duplication)/del (false deletion) for a non-carrier."""
    if carrier:
        return "tp" if cn >= DUP else "fn"
    return "fp" if cn >= DUP else "del" if cn <= DEL else "tn"


def write_tsv(name, header, rows):
    TABLES.mkdir(parents=True, exist_ok=True)
    with open(TABLES / f"{name}.tsv", "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.writelines("\t".join(str(x) for x in r) + "\n" for r in rows)
    print("wrote", TABLES / f"{name}.tsv")
