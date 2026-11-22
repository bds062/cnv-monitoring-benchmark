"""CNVkit duplications against read depth (Additional file 1): every filtered H. zea CNVkit duplication, its
boundaries, the read-depth copy number of the moth inside it (Delly2's `delly cnv` windows of 10 kb mappable sequence, mean +- SE over
the windows) and that of CNVkit's reference (the crosstrio-3 moths of the moth's trio), so that each call can be
checked by eye. Tables: cnvkit_dups_depth.tsv, cnvkit_dups_depth.tex (the supplementary table).
    python analysis/cnvkit_depth.py"""
import bisect
import csv
import math

import numpy as np

from cnvbench import calls as C, locus as L, mir as M
from cnvbench.paths import REPEATS, REPO, TABLES

SP, SV = "zea", "DUP"
ROLES = (("resistant", "R"), ("susceptible", "S"), ("child", "C"))
_win = {}


def windows(s):
    if s not in _win:
        d = {}
        for c, a, b, cn in L.depth_windows(s, SP):
            d.setdefault(c, []).append((a, b, cn))
        _win[s] = {c: tuple(np.array(x) for x in zip(*v)) for c, v in d.items()}
    return _win[s]


def depth(s, c):
    """(mean, SE, windows) of the moth's depth copy number over the call, windows weighted by overlap; SE is nan
    for a call inside one window."""
    st, en, cn = windows(s)[c.chrom]
    i, j = bisect.bisect_right(en, c.start), bisect.bisect_left(st, c.end)
    ov = np.minimum(en[i:j], c.end) - np.maximum(st[i:j], c.start)
    v, ov = cn[i:j][ov > 0], ov[ov > 0]
    se = float(np.std(v, ddof=1) / math.sqrt(len(v))) if len(v) > 1 else math.nan
    return float(np.sum(ov * v) / np.sum(ov)), se, len(v)


def compute():
    trios, rep = C.trios(SP), M.load_repeats(SP, REPEATS[SP])
    filt = M.apply_filters(C.headline(SP)["CNVkit"], trios, rep)
    ref = {}
    for r in csv.reader(open(REPO / "data" / "samples" / "backgrounds.tsv"), delimiter="\t"):
        if r and r[0] == SP and r[1] == "crosstrio3":
            ref.setdefault(r[2], []).append(r[3])
    rows = []
    for t in trios:
        child = [c for c in filt[t.child] if c.svtype == SV]
        cons = M.matched(child, [c for p in t.parents for c in filt[p] if c.svtype == SV])
        for (role, tag), s in zip(ROLES, (t.parents[0], t.parents[1], t.child)):
            for c in sorted((c for c in filt[s] if c.svtype == SV), key=lambda c: (int(c.chrom), c.start)):
                d, se, nw = depth(s, c)
                rd = [depth(x, c) for x in ref[t.name]]
                rows.append({"trio": t.name, "sample": s, "role": tag, "chrom": c.chrom, "start": c.start,
                             "end": c.end, "length_kb": round(M.length(c) / 1e3, 1), "windows": nw,
                             "depth_cn": round(d, 2), "depth_se": round(se, 2),
                             "reference_cn": round(float(np.mean([x[0] for x in rd])), 2),
                             "reference_se": round(float(np.std([x[0] for x in rd], ddof=1) / math.sqrt(len(rd))), 2),
                             "parent_call": ("yes" if cons[child.index(c)] else "no") if role == "child" else "",
                             "tryp77": c.chrom == L.CHROM and c.start < L.AMPLICON[1] and c.end > L.AMPLICON[0]})
    return rows


def pm(x, e):
    return f"{x:.2f}" if math.isnan(e) else f"{x:.2f} $\\pm$ {e:.2f}"


def num(x):
    return f"{x:,}".replace(",", "{,}")


def latex(rows):
    head = ("Trio & Moth & Chr. & Start & End & kb & Win. & Moth CN & Reference CN & Parent call \\\\\n")
    out = ["\\begingroup\\color{red}\\footnotesize\\setlength{\\tabcolsep}{4pt}",
           "\\begin{longtable}{@{}llrrrrrccc@{}}",
           "\\caption{\\bds{\\textbf{Every CNVkit duplication against read depth in the \\textit{H. zea} trios.} "
           "Filtered calls (repeat filter and blacklist, crosstrio-3 reference), "
           f"{len(rows)} in all. Moth: R, resistant parent; S, susceptible parent; C, child. Start, end: call "
           "boundaries. Win.: Delly2 \\texttt{delly cnv} read-depth windows overlapping the call (BWA-MEM "
           "alignments); each window holds 10\\,kb of mappable sequence, so windows are longer in poorly mappable "
           "sequence. Moth CN: the moth's read-depth copy number inside the call, mean of the windows "
           "weighted by overlap $\\pm$ standard error over the windows (no standard error for a call inside one "
           "window). Reference CN: the same for CNVkit's reference, mean $\\pm$ standard error over its three "
           "crosstrio-3 moths. A duplication with a moth CN near 2 over a reference CN near 1 is a deletion "
           "carried by the reference. Parent call: for a child, whether a parent has a matching call "
           "(OC $\\geq 0.5$; no, MIR-inconsistent). $^{\\dagger}$Overlaps the \\textit{Tryp77} ddPCR amplicon. "
           "Every call: \\texttt{cnvkit\\_dups\\_depth.tsv} in the analysis repository.}}",
           "\\label{tab:cnvkitdups}\\\\", "\\toprule", head + "\\midrule", "\\endfirsthead",
           "\\multicolumn{10}{l}{\\textit{(continued)}}\\\\", "\\toprule", head + "\\midrule", "\\endhead",
           "\\bottomrule", "\\endlastfoot"]
    prev = None
    for r in rows:
        if prev and r["trio"] != prev:
            out.append("\\midrule")
        trio = r["trio"].replace("Trio", "") if r["trio"] != prev else ""
        prev = r["trio"]
        dag = "$^{\\dagger}$" if r["tryp77"] else ""
        out.append(f"{trio} & {r['role']} & {r['chrom']}{dag} & {num(r['start'])} & {num(r['end'])} & "
                   f"{r['length_kb']:.1f} & {r['windows']} & {pm(r['depth_cn'], r['depth_se'])} & "
                   f"{pm(r['reference_cn'], r['reference_se'])} & {r['parent_call']} \\\\")
    out += ["\\end{longtable}", "\\endgroup", ""]
    return "\n".join(out)


if __name__ == "__main__":
    rows = compute()
    TABLES.mkdir(parents=True, exist_ok=True)
    with open(TABLES / "cnvkit_dups_depth.tsv", "w") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    open(TABLES / "cnvkit_dups_depth.tex", "w").write(latex(rows))
    print(len(rows), "duplications;", sum(r["windows"] < 2 for r in rows), "inside one window")
