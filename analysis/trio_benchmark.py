"""Trio benchmark: every Mendelian-inconsistency result for both species, computed once from the caller
outputs and cached in results/tables/trio_benchmark.json, the summary tables (TSV) and the MIR figures
(Fig. 2; Additional file 1, Figs. S1-S7, S16 and S20-S24; Tables S2-S4, S7, S8 and S18).
    python analysis/trio_benchmark.py [--recompute]"""
import json
import math
import sys

import matplotlib.lines
import matplotlib.ticker
import matplotlib.pyplot as plt
import matplotlib.transforms
import numpy as np

from cnvbench import calls as C, figstyle, locus as L, mir as M, style
from cnvbench.paths import REF, REPEATS, REPO, TABLES
from cnvbench.style import COLORS, INK, MM, MUTED, ORDER, tint

CACHE = TABLES / "trio_benchmark.json"
SE_CACHE = TABLES / "trio_benchmark_se.json"           # per-trio counts behind the pooled values
INDIV = ["Delly2", "smoove"]
CAP = 1_000_000                                   # parental calls > 1 Mb (Additional file 1, Section S3.3)
SVL = {"DUP": "duplications", "DEL": "deletions"}
RECIP = M.Match("recip", 0.5)
# one setting changed from the headline at a time: (key, label, filter settings, match rule)
SETTINGS = [("headline", "Headline", {}, M.OC),
            ("no_filters", "No filters", {"rep_frac": None, "bl_rule": None}, M.OC),
            ("bl_none", "Blacklist none", {"bl_rule": None}, M.OC),
            ("bl_parents4", "Blacklist ≥4/8 parents", {"bl_rule": ("parents", 4)}, M.OC),
            ("bl_parents5", "Blacklist ≥5/8 parents", {"bl_rule": ("parents", 5)}, M.OC),
            ("bl_parents8", "Blacklist 8/8 parents", {"bl_rule": ("parents", 8)}, M.OC),
            ("bl_trios4", "Blacklist all 4 trios", {"bl_rule": ("trios", 4)}, M.OC),
            ("rep_none", "Repeat filter none", {"rep_frac": None}, M.OC),
            ("rep_70", "Repeat filter 70%", {"rep_frac": 0.7}, M.OC),
            ("rep_25", "Repeat filter 25%", {"rep_frac": 0.25}, M.OC),
            ("oc_80", "OC ≥ 0.8", {}, M.Match("oc", 0.8)),
            ("oc_ratio", "OC + size ratio ≥ 0.1", {}, M.Match("oc", 0.5, 0.1)),
            ("recip", "Reciprocal ≥ 0.5", {}, RECIP)]
POPDESIGN = {"Delly2": ["pooled12", "susceptible13", "neutral18", "neutral9", "crosstrio9"],
             "smoove": ["pooled12", "susceptible13", "neutral18", "neutral9", "crosstrio9"],
             "GATK-gCNV": ["pooled12", "susceptible13", "neutral18", "neutral9", "crosstrio9"],
             "CNVcaller": ["pooled12", "susceptible13", "neutral18", "neutral9", "crosstrio9"],
             "CNVkit": ["pooled12", "susceptible13", "neutral18", "crosstrio3"]}
BG_LABEL = {"pooled12": "pooled-12", "susceptible13": "susceptible-13", "neutral18": "balanced-18",
            "neutral9": "balanced-9", "crosstrio9": "crosstrio-9", "crosstrio3": "crosstrio-3"}


# ─── computation ───────────────────────────────────────────────────────────────────────────────────
def species_results(sp):
    trios = C.trios(sp)
    parents = [p for t in trios for p in t.parents]
    rep = M.load_repeats(sp, REPEATS[sp])
    gsize = M.genome_size(sp, f"{REF[sp]}.fai")
    raw = C.headline(sp)
    out = {"headline": {}, "robust": {}, "stages": {}, "big_calls": []}
    for m in ORDER:
        print(C.SPECIES[sp], m, flush=True)
        filt = {k: M.apply_filters(raw[m], trios, rep, **f) for k, _, f, _ in SETTINGS if f or k == "headline"}
        head, rep_only, bl_only = filt["headline"], filt["bl_none"], filt["rep_none"]
        h, rb, st = {}, {}, {}
        for sv in M.SVTYPES:
            n_all = lambda calls: sum(c.svtype == sv for cs in calls.values() for c in cs)
            h[sv] = {"raw": M.summary(M.per_trio(raw[m], trios, sv)),
                     "filtered": M.summary(M.per_trio(head, trios, sv)),
                     "filtered_recip": M.mir(head, trios, sv, rule=RECIP),
                     "kept": {k: 100 * n_all(v) / n_all(raw[m]) for k, v in
                              (("repeat", rep_only), ("blacklist", bl_only), ("both", head))},
                     "bins": {lab: M.per_trio(head, trios, sv, lo=lo, hi=hi) for lo, hi, lab in M.BINS},
                     "child_len": [M.length(c) for t in trios for c in head[t.child] if c.svtype == sv]}
            rb[sv] = {}
            for k, _, f, rule in SETTINGS:
                pt = M.per_trio(filt.get(k, head), trios, sv, rule=rule)
                rb[sv][k] = [M.pct(sum(x[0] for x in pt), sum(x[1] for x in pt)), sum(x[1] for x in pt)]
            st[sv] = {}
            for stage, calls, cap in (("raw", raw[m], None), ("repeat", rep_only, None),
                                      ("blacklist", head, None), ("cap1mb", head, CAP)):
                typed = {s: [c for c in cs if c.svtype == sv and (cap is None or s not in parents or M.length(c) <= cap)]
                         for s, cs in calls.items()}
                pt = M.per_trio(typed, trios, sv)
                st[sv][stage] = {"mir": M.pct(sum(x[0] for x in pt), sum(x[1] for x in pt)),
                                 "n_child": sum(x[1] for x in pt),
                                 "n_big_parent": sum(M.length(c) > CAP for p in parents for c in typed[p]),
                                 "cover": [M.cover(typed[p], gsize) if typed[p] else 0.0 for p in parents]}
            if sp == "zea" and m in INDIV:
                out["big_calls"] += big_calls(m, head, trios, sv)
        out["headline"][m], out["robust"][m], out["stages"][m] = h, rb, st
    return out


def big_calls(m, head, trios, sv):
    """Filtered parental calls > 1 Mb: depth copy number of the parent inside the call (`delly cnv`) and the
    child calls whose only parental match it is."""
    out = []
    for t in trios:
        child = [c for c in head[t.child] if c.svtype == sv]
        small = [c for p in t.parents for c in head[p] if c.svtype == sv and M.length(c) <= CAP]
        unmatched = [c for c, ok in zip(child, M.matched(child, small)) if not ok]
        for p in t.parents:
            for c in head[p]:
                if c.svtype == sv and M.length(c) > CAP:
                    only = [k for k, ok in zip(unmatched, M.matched(unmatched, [c])) if ok]
                    out.append({"caller": m, "trio": t.name, "parent": p, "resistant": p == t.parents[0],
                                "chrom": c.chrom, "start": c.start, "end": c.end, "svtype": sv,
                                "depth_cn": L.caller_cn("Delly2", "none", p, c.start, c.end, c.chrom),
                                "only_child": [(k.start, k.end) for k in only]})
    return out


def popdesign():
    """H. zea callers under each background: filtered MIR per trio (blacklist built
    within each background's own calls)."""
    trios = C.trios("zea")
    rep = M.load_repeats("zea", REPEATS["zea"])
    out = {}
    for m, bgs in POPDESIGN.items():
        out[m] = {}
        for bg in bgs:
            print("H. zea", m, bg, flush=True)
            raw = C.load_set("zea", {bg: (m, bg)})[bg]
            filt = M.apply_filters(raw, trios, rep)
            out[m][bg] = {sv: M.per_trio(filt, trios, sv) for sv in M.SVTYPES}
    return out


def cnvkit_protocols():
    """CNVkit raw and filtered MIR under the earlier protocol (crosstrio3_auto: CNVkit's automatic WGS bins,
    4,589 bp in H. zea) and the matched headline protocol (crosstrio3), both species."""
    out = {}
    for sp in ("zea", "human"):
        trios, rep = C.trios(sp), M.load_repeats(sp, REPEATS[sp])
        out[sp] = {}
        for bg in ("crosstrio3_auto", "crosstrio3"):
            raw = C.load_set(sp, {bg: ("CNVkit", bg)})[bg]
            filt = M.apply_filters(raw, trios, rep)
            out[sp][bg] = {sv: {cond: M.summary(M.per_trio(calls, trios, sv)) for cond, calls in
                                (("raw", raw), ("filtered", filt))} for sv in M.SVTYPES}
    return out


def se_extras(sp):
    """Per-trio counts behind the pooled values of the figures and the paired caller comparisons:
    [(inconsistent, child calls)] per trio for each filter stage, robustness setting and the reciprocal rule,
    and [(calls kept, all calls)] per trio for each filter (the three moths of a trio pooled)."""
    trios = C.trios(sp)
    parents = [p for t in trios for p in t.parents]
    rep = M.load_repeats(sp, REPEATS[sp])
    raw = C.headline(sp)
    out = {}
    for m in ORDER:
        print("SE", C.SPECIES[sp], m, flush=True)
        filt = {k: M.apply_filters(raw[m], trios, rep, **f) for k, _, f, _ in SETTINGS if f or k == "headline"}
        head = filt["headline"]
        o = {}
        for sv in M.SVTYPES:
            st = {}
            for stage, calls, cap in (("raw", raw[m], None), ("repeat", filt["bl_none"], None),
                                      ("blacklist", head, None), ("cap1mb", head, CAP)):
                typed = {s: [c for c in cs if c.svtype == sv and (cap is None or s not in parents or M.length(c) <= cap)]
                         for s, cs in calls.items()}
                st[stage] = M.per_trio(typed, trios, sv)
            n = lambda calls, t: sum(c.svtype == sv for s in (t.child, *t.parents) for c in calls[s])
            o[sv] = {"stages": st,
                     "robust": {k: M.per_trio(filt.get(k, head), trios, sv, rule=rule) for k, _, f, rule in SETTINGS},
                     "recip": M.per_trio(head, trios, sv, rule=RECIP),
                     "kept": {k: [(n(v, t), n(raw[m], t)) for t in trios] for k, v in
                              (("repeat", filt["bl_none"]), ("blacklist", filt["rep_none"]), ("both", head))}}
        out[m] = o
    return out


def load_se(R):
    """Attach the per-trio counts (cached in trio_benchmark_se.json) and check them against the main cache."""
    if SE_CACHE.exists() and "--recompute-se" not in sys.argv:
        S = json.load(open(SE_CACHE))
    else:
        S = {sp: se_extras(sp) for sp in ("zea", "human")}
        json.dump(S, open(SE_CACHE, "w"))
    for sp in ("zea", "human"):
        for m in ORDER:
            for sv in M.SVTYPES:
                for stage, pt in S[sp][m][sv]["stages"].items():
                    a, b = pooled(pt), R[sp]["stages"][m][sv][stage]["mir"]
                    assert (a is None and b is None) or abs(a - b) < 1e-9, (sp, m, sv, stage, a, b)
    R["se"] = S
    return R


def trio_vals(pt):
    """Per-trio MIR (%) from [(inconsistent, child calls)]; trios without child calls are left out."""
    return [100 * k / n for k, n in pt if n]


def trio_dots(ax, pos, vals, color, horizontal=False, s=5, hollow=False, marker="o", z=2.5, lw=0.5):
    """The per-trio (or per-parent) values behind a pooled value, as small dots at one x (or y) position."""
    v = [x for x in vals if x is not None and np.isfinite(x)]
    xy = (v, [pos] * len(v)) if horizontal else ([pos] * len(v), v)
    kw = dict(facecolor="white", edgecolor=color, lw=lw) if hollow else dict(color=color, lw=0, alpha=0.7)
    ax.scatter(*xy, s=s, marker=marker, zorder=z, **kw)


def compute():
    res = {sp: species_results(sp) for sp in ("zea", "human")}
    res["popdesign"] = popdesign()
    res["cnvkit_protocol"] = cnvkit_protocols()
    res["tryp_depth_cn"] = {t.parents[0]: L.caller_cn("Delly2", "none", t.parents[0], *L.TRYP) for t in C.trios("zea")}
    TABLES.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(CACHE, "w"))
    return res


# ─── tables ────────────────────────────────────────────────────────────────────────────────────────
def f1(x):
    return "" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.1f}"


def write_tables(R):
    with open(TABLES / "mir_headline.tsv", "w") as fh:          # Table S3 (per-trio MIR)
        fh.write("species\tsvtype\tcaller\tcalls\tn\tinconsistent\tmir\t"
                 "per_trio\tmir_reciprocal\tkept_repeat\tkept_blacklist\tkept_both\n")
        for sp in ("zea", "human"):
            for sv in M.SVTYPES:
                for m in ORDER:
                    h = R[sp]["headline"][m][sv]
                    for cond in ("raw", "filtered"):
                        s = h[cond]
                        extra = [f1(h["filtered_recip"])] + [f1(h["kept"][k]) for k in ("repeat", "blacklist", "both")] \
                            if cond == "filtered" else [""] * 4
                        fh.write("\t".join([C.SPECIES[sp], sv, m, cond, str(s["n"]), str(s["inc"]), f1(s["mir"]),
                                            ", ".join(map(f1, s["per_trio"])), *extra]) + "\n")
    with open(TABLES / "mir_paired.tsv", "w") as fh:            # Table S4
        # every pair of callers, filtered main-analysis MIR, compared within trios (cnvbench.mir.paired)
        fh.write("species\tsvtype\tcaller_a\tcaller_b\tmir_a\tmir_b\tb_above_a_trios\ttrios\todds_ratio_b_vs_a\tp\n")
        for sp in ("zea", "human"):
            for sv in M.SVTYPES:
                for i, a in enumerate(ORDER):
                    for b in ORDER[i + 1:]:
                        pa, pb = (R["se"][sp][m][sv]["robust"]["headline"] for m in (a, b))
                        r = M.paired(pa, pb)
                        fh.write(f"{C.SPECIES[sp]}\t{sv}\t{a}\t{b}\t{f1(pooled(pa))}\t{f1(pooled(pb))}\t{r['above']}\t"
                                 f"{r['trios']}\t{r['odds_ratio']:.2f}\t{r['p']:.3g}\n")
    with open(TABLES / "mir_robustness.tsv", "w") as fh:        # Tables S18 (duplications) and S2 (deletions)
        fh.write("svtype\tsetting\t" + "\t".join(f"{C.SPECIES[sp]} {m}" for sp in ("zea", "human") for m in ORDER) + "\n")
        for sv in M.SVTYPES:
            for k, lab, _, _ in SETTINGS:
                fh.write(f"{sv}\t{lab}\t" + "\t".join(f1(R[sp]["robust"][m][sv][k][0]) for sp in ("zea", "human")
                                                      for m in ORDER) + "\n")
    with open(TABLES / "mir_size_bins.tsv", "w") as fh:
        fh.write("species\tsvtype\tcaller\tsize_bin\tn\tmir\n")
        for sp in ("zea", "human"):
            for sv in M.SVTYPES:
                for m in ORDER:
                    for _, _, lab in M.BINS:
                        pt = R[sp]["headline"][m][sv]["bins"][lab]
                        n = sum(x[1] for x in pt)
                        fh.write(f"{C.SPECIES[sp]}\t{sv}\t{m}\t{lab}\t{n}\t{f1(M.pct(sum(x[0] for x in pt), n))}\n")
    with open(TABLES / "mir_popdesign.tsv", "w") as fh:         # Figs. S4 and S22
        names = [t.name for t in C.trios("zea")]
        fh.write("caller\tbackground\tsvtype\t" + "\t".join(f"mir_{t}\tn_{t}" for t in names) + "\tmean_mir\tpooled_mir\n")
        for m, d in R["popdesign"].items():
            for bg, v in d.items():
                for sv, pt in v.items():
                    per = [M.pct(*x) for x in pt]
                    fh.write(f"{m}\t{BG_LABEL[bg]}\t{sv}\t" + "\t".join(f"{f1(p)}\t{n}" for p, (_, n) in zip(per, pt))
                             + f"\t{f1(np.nanmean([np.nan if p is None else p for p in per]))}"
                             + f"\t{f1(M.pct(sum(x[0] for x in pt), sum(x[1] for x in pt)))}\n")
    with open(TABLES / "mir_cnvkit_protocol.tsv", "w") as fh:   # Table S7
        fh.write("species\tprotocol\tsvtype\tcalls\tmir\tinconsistent\tn\n")
        for sp, d in R["cnvkit_protocol"].items():
            for bg, lab in (("crosstrio3_auto", "earlier"), ("crosstrio3", "matched")):
                for sv in M.SVTYPES:
                    for cond, s in d[bg][sv].items():
                        fh.write(f"{C.SPECIES[sp]}\t{lab}\t{sv}\t{cond}\t{f1(s['mir'])}\t{s['inc']}\t{s['n']}\n")
    with open(TABLES / "mir_oversized.tsv", "w") as fh:         # Table S8
        fh.write("species\tcaller\tsvtype\tmir_headline\tmir_parent_le_1mb\tparent_calls_gt_1mb\n")
        for sp in ("zea", "human"):
            for m in ORDER:
                for sv in M.SVTYPES:
                    s = R[sp]["stages"][m][sv]
                    fh.write(f"{C.SPECIES[sp]}\t{m}\t{sv}\t{f1(s['blacklist']['mir'])}\t{f1(s['cap1mb']['mir'])}\t"
                             f"{s['blacklist']['n_big_parent']}\n")
    with open(TABLES / "filter_stages.tsv", "w") as fh:         # Fig. 2 and Fig. S20
        fh.write("species\tcaller\tsvtype\tstage\tmir\tn_child\tparent_cover_median\tparent_cover\n")
        for sp in ("zea", "human"):
            for m in ORDER:
                for sv in M.SVTYPES:
                    for stage, s in R[sp]["stages"][m][sv].items():
                        fh.write(f"{C.SPECIES[sp]}\t{m}\t{sv}\t{stage}\t{f1(s['mir'])}\t{s['n_child']}\t"
                                 f"{np.median(s['cover']):.2f}\t" + ", ".join(f"{x:.2f}" for x in s["cover"]) + "\n")


# ─── figures ───────────────────────────────────────────────────────────────────────────────────────
def pooled(pt):
    return M.pct(sum(x[0] for x in pt), sum(x[1] for x in pt))


def nz(x):
    return np.nan if x is None else x


def eb(ax, x, y, err, color, horizontal=False, lw=0.8):
    """Error bars (no markers): err is ±1 standard error, or [[below], [above]] lengths of an interval."""
    kw = dict(fmt="none", ecolor=color, elinewidth=lw, capsize=1.6, capthick=lw, zorder=2.5)
    ax.errorbar(x, y, xerr=err, **kw) if horizontal else ax.errorbar(x, y, yerr=err, **kw)


def binom_se(k, n):
    return 100 * math.sqrt((k / n) * (1 - k / n) / n) if n else np.nan


def group_split_x(ax, y=-0.4):
    ax.axvline(len(INDIV) - 0.5, color=MUTED, lw=0.7, ls=":", zorder=1)
    tr = matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes)
    for lab, x in (("Breakpoint", (len(INDIV) - 1) / 2), ("Read depth", (len(INDIV) + len(ORDER) - 1) / 2)):
        ax.text(x, y, lab, transform=tr, ha="center", va="top", fontsize=6.3, color=MUTED)


def group_split_y(ax, ys, labels=True):
    ax.axhline((ys[len(INDIV) - 1] + ys[len(INDIV)]) / 2, color=MUTED, lw=0.7, ls=":", zorder=1)
    if labels:
        tr = matplotlib.transforms.blended_transform_factory(ax.transAxes, ax.transData)
        for lab, rows in (("Breakpoint", ys[:len(INDIV)]), ("Read depth", ys[len(INDIV):])):
            ax.text(1.02, np.mean(rows), lab, transform=tr, rotation=270, ha="left", va="center", fontsize=6.3,
                    color=MUTED)


def fig_metric_filter(R, sv):
    Z = R["zea"]["headline"]
    greys = ["#C9C9C9", "#8C8C8C", "#3F3F3F"]
    fig, axes = plt.subplots(1, 3, figsize=(170 * MM, 72 * MM), gridspec_kw={"width_ratios": [1.15, 1.15, 1]})
    x = np.arange(len(ORDER))
    ax = axes[0]
    SE = R["se"]["zea"]
    for i, (vals, per, lab) in enumerate(((lambda m: Z[m][sv]["filtered_recip"], lambda m: SE[m][sv]["recip"],
                                           "50% reciprocal overlap"),
                                          (lambda m: Z[m][sv]["filtered"]["mir"],
                                           lambda m: SE[m][sv]["stages"]["blacklist"], "Overlap coef. ≥ 0.5"))):
        xi = x + (i - 0.5) * 0.36
        ax.bar(xi, [nz(vals(m)) for m in ORDER], 0.33, color=greys[2 * i], label=lab, zorder=2)
        for xx, m in zip(xi, ORDER):
            trio_dots(ax, xx, trio_vals(per(m)), INK, s=6, hollow=True, lw=0.6, z=3)
    ax.set(ylim=(0, 100), ylabel="Filtered MIR (%)")
    ax.set_title("Match rule", loc="left")
    ax = axes[1]
    for i, (k, lab) in enumerate((("repeat", "Repeat filter"), ("blacklist", "Blacklist"), ("both", "Both"))):
        xi = x + (i - 1) * 0.26
        ax.bar(xi, [Z[m][sv]["kept"][k] for m in ORDER], 0.24, color=greys[i], label=lab, zorder=2)
        for xx, m in zip(xi, ORDER):
            trio_dots(ax, xx, trio_vals(SE[m][sv]["kept"][k]), INK, s=6, hollow=True, lw=0.6, z=3)
    ax.set(ylim=(0, 130), yticks=[0, 20, 40, 60, 80, 100], ylabel="Calls kept (%)")
    ax.set_title("Calls kept by each filter", loc="left")
    for ax, k in zip(axes[:2], "ab"):
        ax.set_xticks(x, ORDER, rotation=30, ha="right")
        ax.spines["left"].set_bounds(0, 100)
        ax.grid(axis="y", zorder=0)
        ax.legend(frameon=False, loc="upper left")
        group_split_x(ax)
        style.label(ax, k, x=-0.14)
    ax = axes[2]
    ys = np.arange(len(ORDER))[::-1]
    for y, m in zip(ys, ORDER):
        r0, r1 = Z[m][sv]["raw"]["mir"], Z[m][sv]["filtered"]["mir"]
        ax.plot([r0, r1], [y, y], color="#BBBBBB", lw=1, zorder=1)
        trio_dots(ax, y + 0.22, trio_vals(SE[m][sv]["stages"]["raw"]), COLORS[m], horizontal=True, s=5, hollow=True)
        trio_dots(ax, y - 0.22, trio_vals(SE[m][sv]["stages"]["blacklist"]), COLORS[m], horizontal=True, s=5)
        ax.scatter([r0], [y], s=20, facecolor="white", edgecolor=COLORS[m], lw=1.1, zorder=2)
        ax.scatter([r1], [y], s=20, color=COLORS[m], zorder=3)
    ax.set_yticks(ys, ORDER)
    ax.set_ylim(-0.6, len(ORDER) + 0.9)
    ax.set_xlabel("MIR, overlap coef. ≥ 0.5 (%)")
    ax.scatter([], [], s=20, facecolor="white", edgecolor="#666666", lw=1.1, label="Unfiltered")
    ax.scatter([], [], s=20, color="#666666", label="Filtered")
    ax.legend(frameon=False, loc="upper right", ncol=2, handletextpad=0.2, columnspacing=0.8)
    ax.set_title("Effect of filtering", loc="left")
    ax.grid(axis="x", zorder=0)
    group_split_y(ax, ys, labels=False)
    style.label(ax, "c", x=-0.14)
    fig.tight_layout(w_pad=1.2)
    style.save(fig, "fig_metric_filter_" + sv.lower())


def gene_lengths():
    rows = [l.split("\t") for l in open(REPO / "data" / "ogs2_gene_lengths.tsv") if not l.startswith("#")]
    return np.array([int(r[1]) for r in rows[1:]])


def fig_size_bins(R, sv):
    Z = R["zea"]["headline"]
    fig, axes = plt.subplots(1, 2, figsize=(170 * MM, 58 * MM), gridspec_kw={"width_ratios": [1, 1.3]})
    g = gene_lengths()
    ax = axes[0]
    edges = [2, np.log10(2_000), 4, np.log10(50_000), 6]
    cols = ["#C9D4E0", "#93A8BF", "#5E7C9C", "#2F4B6B"]
    labs = ["0.1–2 kb", "2–10", "10–50", "≥50"]
    counts, e = np.histogram(np.log10(g), bins=np.linspace(1, 6, 81))
    mids = (e[:-1] + e[1:]) / 2
    for k in range(4):
        sel = (mids >= edges[k]) & (mids < edges[k + 1])
        ax.bar(mids[sel], counts[sel], width=np.diff(e)[0] * 0.92, color=cols[k], zorder=2)
    ax.bar(mids[mids < edges[0]], counts[mids < edges[0]], width=np.diff(e)[0] * 0.92, color="#E4E4E4", zorder=2)
    top = counts.max() * 1.4
    shares = [((g >= lo) & (g < hi)).mean() for lo, hi, _ in M.BINS]
    for k in range(4):
        x0, x1 = edges[k], edges[k + 1]
        ax.plot([x0 + 0.02, x1 - 0.02], [top * 0.86] * 2, color=cols[k], lw=2.2, solid_capstyle="butt")
        ax.text((x0 + x1) / 2, top * 0.89, labs[k], ha="center", va="bottom", fontsize=6.2, color="#333333")
        ax.text((x0 + x1) / 2, top * 0.83, f"{shares[k] * 100:.0f}%", ha="center", va="top", fontsize=6.2,
                color="#333333")
    med = np.median(g)
    ax.plot([np.log10(med)] * 2, [0, counts.max() * 1.02], color="#C62828", lw=0.9, zorder=3)
    ax.text(np.log10(med) + 0.06, counts.max() * 0.98, f"median {med / 1000:.1f} kb", color="#C62828", fontsize=6.2,
            va="top")
    ticks = [100, 1_000, 10_000, 100_000, 1_000_000]
    names = ["100 bp", "1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb"]
    ymax = np.floor(counts.max() / 200) * 200
    ax.set_xticks(np.log10(ticks), names[:5])
    ax.set(xlim=(1.6, 6.05), ylim=(0, top), xlabel=f"Gene length (lifted OGS2, n = {len(g):,})", ylabel="Genes",
           yticks=np.arange(0, ymax + 1, 200))
    ax.spines["left"].set_bounds(0, ymax)
    ax.grid(axis="y", zorder=0)
    ax.set_title(r"$\it{H.\ zea}$ gene lengths and MIR size bins", loc="left")
    style.label(ax, "a", x=-0.14)
    ax = axes[1]
    parts = ax.violinplot([np.log10(Z[m][sv]["child_len"]) for m in ORDER], positions=range(len(ORDER)), vert=False,
                          showmedians=True, widths=0.8)
    for body, m in zip(parts["bodies"], ORDER):
        body.set(facecolor=COLORS[m], edgecolor="white", alpha=0.85)
    for k in ("cmedians", "cmins", "cmaxes", "cbars"):
        parts[k].set(color="#333333", linewidth=0.6)
    for lo, _, _ in M.BINS:
        ax.axvline(np.log10(lo), color="#9A9A9A", lw=0.6, ls=(0, (3, 2)), zorder=0)
    ax.set_yticks(range(len(ORDER)), [f"{m} (n={len(Z[m][sv]['child_len']):,})" for m in ORDER])
    ax.invert_yaxis()
    group_split_y(ax, list(range(len(ORDER))))
    ax.set_xticks(np.log10(ticks + [10_000_000]), names)
    ax.set_xlabel(f"Child call length (filtered {SVL[sv]})")
    ax.set_title(SVL[sv].capitalize() + r" by caller, $\it{H.\ zea}$", loc="left")
    style.label(ax, "b", x=-0.14)
    fig.tight_layout(w_pad=1.5)
    style.save(fig, "fig_size_bins_" + sv.lower())


def msize(n):
    return 9 + 200 * np.asarray(n, dtype=float) / 5000           # marker area ~ child calls


def fig_mir_by_size(R, sv):
    fig, axes = plt.subplots(1, 2, figsize=(170 * MM, 62 * MM), sharey=True)
    labs = [b[2] for b in M.BINS]
    xs = np.arange(len(labs))
    for ax, sp, k in zip(axes, ("zea", "human"), "ab"):
        for m in ORDER:
            pts = [R[sp]["headline"][m][sv]["bins"][b] for b in labs]
            n = np.array([sum(x[1] for x in pt) for pt in pts])
            y = np.array([nz(pooled(pt)) for pt in pts])
            ok = n > 0                                            # bins without calls are left out
            ax.plot(xs[ok], y[ok], color=COLORS[m], lw=1.6, label=m, zorder=3, ls=(0, (3, 1.5)) if m in INDIV else "-")
            ax.scatter(xs[ok], y[ok], s=msize(n[ok]), color=COLORS[m], edgecolor="white", lw=0.5, zorder=4, alpha=0.9)
        ax.set_xticks(xs, ["0.1–2 kb", "2–10 kb", "10–50 kb", "≥50 kb"])
        ax.set_xlim(-0.3, len(labs) - 0.2)
        ax.set_title(C.SPECIES[sp], loc="left", style="italic" if sp == "zea" else "normal")
        ax.set_xlabel("Child call length")
        ax.grid(axis="y", zorder=0)
        ax.text(-0.02, 1.04, k, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right")
    axes[0].set_ylabel(f"Filtered MIR, {SVL[sv]} (%)")
    sz = [axes[1].scatter([], [], s=msize([n]), color="#9A9A9A", edgecolor="white", lw=0.5) for n in (10, 100, 1000, 5000)]
    axes[1].legend(sz, ["10", "100", "1,000", "5,000"], title="Child calls (marker area)", title_fontsize=6.3,
                   frameon=False, loc="upper right", ncol=4, handletextpad=0.1, columnspacing=0.9, borderaxespad=0.2)
    hl = dict(zip(*axes[0].get_legend_handles_labels()[::-1]))
    fig.tight_layout(w_pad=1.5, rect=(0, 0.1, 1, 1))
    for grp, title, x in ((INDIV, "Breakpoint (dashed)", 0.28),
                          ([m for m in ORDER if m not in INDIV], "Read depth (solid)", 0.68)):
        fig.legend([hl[m] for m in grp], grp, frameon=False, loc="lower center", ncol=len(grp), bbox_to_anchor=(x, 0.0),
                   title=title, title_fontsize=6.5)
    style.save(fig, "fig_mir_by_size_" + sv.lower())


def fig_popdesign(R, sv):
    shade = {"susceptible13": 0.25, "neutral18": 0.45, "neutral9": 0.7, "crosstrio9": 1.0, "crosstrio3": 1.0}
    # breakpoint callers on the top row, read-depth callers below
    fig, grid = plt.subplots(2, 3, figsize=(170 * MM, 150 * MM), sharey=True)
    grid[0, 2].axis("off")
    axes = [grid[0, 0], grid[0, 1], grid[1, 0], grid[1, 1], grid[1, 2]]
    names = [t.name for t in C.trios("zea")]
    for c, (ax, m) in enumerate(zip(axes, POPDESIGN)):
        bgs = POPDESIGN[m]
        bw = 0.8 / len(bgs)
        for j, bg in enumerate(bgs):
            pts = R["popdesign"][m][bg][sv]
            per = [nz(M.pct(*x)) for x in pts]
            col = "#BDBDBD" if bg == "pooled12" else tint(COLORS[m], shade[bg])
            mean = np.nanmean(per)
            xb = np.arange(4) + (j - (len(bgs) - 1) / 2) * bw
            ax.bar(xb, per, bw * 0.9, color=col, zorder=2, label=f"{BG_LABEL[bg]}  ({mean:.1f}%)")
            eb(ax, xb, per, [binom_se(*x) for x in pts], "#555555", lw=0.6)
            ax.axhline(mean, color="#8A8A8A" if bg == "pooled12" else col, ls=(0, (1, 1.5)), lw=1.0, zorder=3)
        ax.set_xticks(range(4), names)
        ax.grid(axis="y", zorder=0)
        ax.set(ylim=(0, 100), yticks=[0, 25, 50, 75, 100])
        ax.set_title(m, loc="left")
        if c in (0, 2):
            ax.set_ylabel(f"Filtered MIR, overlap coef.\n$\\geq$ 0.5, {SVL[sv]} (%)")
        else:
            ax.tick_params(axis="y", length=0)
        style.label(ax, "abcde"[c], x=-0.17 if c in (0, 2) else -0.1)
        ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, -0.14), title="Background (mean MIR)",
                  title_fontsize=6.3, alignment="left", handlelength=1.2)
    for ax, lab in ((grid[0, 0], "Breakpoint callers"), (grid[1, 0], "Read-depth callers")):
        ax.text(-0.17, 1.17, lab, transform=ax.transAxes, fontsize=7.5, fontweight="bold", color=MUTED)
    fig.tight_layout(w_pad=1.2, h_pad=4.0)
    style.save(fig, "fig_popdesign_" + sv.lower())


def fig_robustness(R, sv):
    val = lambda sp, k, m: nz(R[sp]["robust"][m][sv][k][0])
    allv = [val(sp, k, m) for sp in ("zea", "human") for k, *_ in SETTINGS for m in ORDER]
    xmax = int(np.ceil(np.nanmax(allv) / 20) * 20) + 5 if np.nanmax(allv) < 80 else 105
    fig, axes = plt.subplots(1, 5, figsize=(170 * MM, 80 * MM), sharey=True)
    ys = np.arange(len(SETTINGS))[::-1]
    for i, (ax, m) in enumerate(zip(axes, ORDER)):
        ax.axhspan(ys[0] - 0.5, ys[0] + 0.5, color="#F1F1F1", lw=0, zorder=0)
        for y, (k, *_) in zip(ys, SETTINGS):
            z, h = val("zea", k, m), val("human", k, m)
            ax.plot([z, h], [y, y], color="#C8C8C8", lw=0.9, zorder=1)
            trio_dots(ax, y + 0.2, trio_vals(R["se"]["zea"][m][sv]["robust"][k]), COLORS[m], horizontal=True, s=3)
            trio_dots(ax, y - 0.2, trio_vals(R["se"]["human"][m][sv]["robust"][k]), COLORS[m], horizontal=True, s=3,
                      hollow=True, marker="s", lw=0.4)
            ax.scatter([h], [y], marker="s", s=14, facecolor="white", edgecolor=COLORS[m], lw=1.0, zorder=3)
            ax.scatter([z], [y], marker="o", s=16, color=COLORS[m], zorder=3)
        ax.axvline(val("zea", "headline", m), color=COLORS[m], lw=0.5, ls="--", alpha=0.6, zorder=0)
        ax.set_xlim(0, xmax)
        ax.set_xticks(np.arange(0, xmax + 1, 20 if xmax <= 80 else 25))
        ax.grid(axis="x", zorder=0)
        ax.set_title(m, loc="left", fontsize=7)
        ax.tick_params(axis="y", length=0)
        if i:
            ax.spines["left"].set_visible(False)
    axes[0].set_yticks(ys, [lab for _, lab, *_ in SETTINGS])
    axes[0].get_yticklabels()[0].set_fontweight("bold")
    axes[2].set_xlabel(f"MIR, {SVL[sv]} (%)")
    axes[-1].scatter([], [], marker="o", s=16, color="#666666", label=r"$\it{H.\ zea}$")
    axes[-1].scatter([], [], marker="s", s=14, facecolor="white", edgecolor="#666666", lw=1.0, label="Human")
    fig.tight_layout(w_pad=0.6, rect=(0, 0.06, 1, 0.93))
    fig.legend(*axes[-1].get_legend_handles_labels(), frameon=False, loc="lower center", ncol=2,
               bbox_to_anchor=(0.6, 0.0))
    p1, p2 = axes[1].get_position(), axes[2].get_position()
    xdiv, top = (p1.x1 + p2.x0) / 2, p1.y1 + 0.075
    fig.add_artist(matplotlib.lines.Line2D([xdiv, xdiv], [p1.y0 - 0.02, top + 0.03], transform=fig.transFigure,
                                           color=MUTED, lw=0.7, ls=":"))
    fig.text((axes[0].get_position().x0 + p1.x1) / 2, top, "Breakpoint", ha="center", fontsize=7.5,
             fontweight="bold", color=MUTED)
    fig.text((p2.x0 + axes[-1].get_position().x1) / 2, top, "Read depth", ha="center", fontsize=7.5,
             fontweight="bold", color=MUTED)
    style.save(fig, "fig_robustness_" + sv.lower())


def fig_species(R, sv):
    fig, ax = plt.subplots(figsize=(85 * MM, 55 * MM))
    ys = np.arange(len(ORDER))[::-1]
    for y, m in zip(ys, ORDER):
        z, h = (R[sp]["headline"][m][sv]["filtered"]["mir"] for sp in ("zea", "human"))
        ax.plot([z, h], [y, y], color="#BBBBBB", lw=1, zorder=1)
        trio_dots(ax, y + 0.2, trio_vals(R["se"]["zea"][m][sv]["stages"]["blacklist"]), COLORS[m], horizontal=True, s=5)
        trio_dots(ax, y - 0.2, trio_vals(R["se"]["human"][m][sv]["stages"]["blacklist"]), COLORS[m], horizontal=True,
                  s=5, hollow=True, marker="s")
        ax.scatter([z], [y], s=22, marker="o", color=COLORS[m], zorder=3)
        ax.scatter([h], [y], s=22, marker="s", facecolor="white", edgecolor=COLORS[m], lw=1.2, zorder=3)
    ax.set_yticks(ys, ORDER)
    group_split_y(ax, ys)
    ax.set_xlabel(f"Filtered MIR, {SVL[sv]}, overlap coef. ≥ 0.5 (%)")
    ax.scatter([], [], marker="o", color="#555555", s=18, label=r"$\it{H.\ zea}$")
    ax.scatter([], [], marker="s", facecolor="white", edgecolor="#555555", s=18, label="Human")
    ax.legend(frameon=False, loc="upper right")
    ax.grid(axis="x", zorder=0)
    fig.tight_layout()
    style.save(fig, "fig_species_" + sv.lower(), max_mm=None)


def fig_filter_stages(R, sv):
    stages = [("raw", "Raw"), ("repeat", "+ repeat\nfilter"), ("blacklist", "+ blacklist"),
              ("cap1mb", "+ parental calls\n> 1 Mb removed")]
    fig, axs = plt.subplots(3, 2, figsize=(170 * MM, 140 * MM), sharex=True, sharey="row")
    fig.subplots_adjust(left=0.08, right=0.99, top=0.955, bottom=0.11, hspace=0.2, wspace=0.08)
    x = np.arange(len(stages))
    for j, sp in enumerate(("zea", "human")):
        ac, an, am = axs[0, j], axs[1, j], axs[2, j]          # cover, child calls, MIR
        for k, m in enumerate(ORDER):
            S, off = R[sp]["stages"][m][sv], (k - 2) * 0.05
            for i, (s, _) in enumerate(stages):
                ac.scatter(np.full(8, i + off), S[s]["cover"], s=5, color=tint(COLORS[m], 0.45), lw=0, zorder=2)
            ac.plot(x + off, [np.median(S[s]["cover"]) for s, _ in stages], color=COLORS[m], lw=1.3, marker="o",
                    ms=3.2, zorder=3, label=m)
            an.plot(x + off, [S[s]["n_child"] for s, _ in stages], color=COLORS[m], lw=1.3, marker="o", ms=3.2, zorder=3)
            am.plot(x + off, [S[s]["mir"] for s, _ in stages], color=COLORS[m], lw=1.3, marker="o", ms=3.2, zorder=3)
        for a in (ac, an, am):
            a.grid(axis="y", zorder=0)
            a.set_axisbelow(True)
        ac.set(ylim=(0, 100), yticks=[0, 25, 50, 75, 100])
        top = max(R[q]["stages"][m][sv][s]["n_child"] for q in ("zea", "human") for m in ORDER for s, _ in stages)
        an.set(yscale="log", ylim=(10, max(1e5, 10 ** np.ceil(np.log10(top)))))
        am.set_ylim(0, max(50, 10 * np.ceil(max(R[q]["stages"][m][sv][s]["mir"] for q in ("zea", "human")
                                                   for m in ORDER for s, _ in stages) / 10)))
        ac.set_title((r"$\it{H.\ zea}$ " if sp == "zea" else "Human ") + SVL[sv], loc="left")
        am.set_xticks(x, [s for _, s in stages])
        style.label(ac, "ab"[j], x=-0.1 if j else -0.13)
        style.label(an, "cd"[j], x=-0.1 if j else -0.13)
        style.label(am, "ef"[j], x=-0.1 if j else -0.13)
    axs[0, 0].set_ylabel("Parent's genome covered\nby its calls (%)")
    axs[1, 0].set_ylabel(f"Child {SVL[sv][:-1]}\ncalls")
    axs[2, 0].set_ylabel("MIR (%)")
    fig.legend(*axs[0, 0].get_legend_handles_labels(), loc="lower center", ncol=5, frameon=False,
               bbox_to_anchor=(0.53, 0.0))
    style.save(fig, "filter_stages_zea_vs_human" + ("" if sv == "DUP" else "_del"))


def fig_filters_and_size(R, sv):
    """Main-text figure (print style, 170 mm): the filter stages (rows 1-3) and MIR by child call length
    (row 4), both species. Stage 'blacklist' is the main analysis (shaded)."""
    F = figstyle
    stages = [("raw", "Raw calls"), ("repeat", "+ Repeat\nfilter"), ("blacklist", "+ Blacklist"),
              ("cap1mb", "+ Ignore\nparent calls\n> 1 Mb")]
    W, H = F.FULL, 194.0
    L_, R_, GAP = 20.5, 1.0, 6.0                                  # mm: left margin, right margin, column gap
    w = (W - L_ - R_ - GAP) / 2
    T_, RH, RG = 6.0, 35.0, 3.0                                   # top margin, height and gap of rows 1-3
    ZH, ZB = 44.0, 17.5                                           # row 4: height and bottom edge
    with F.paper():
        fig = plt.figure(figsize=(W * F.MM, H * F.MM))
        axs = np.empty((4, 2), dtype=object)
        for j in range(2):
            x0 = (L_ + j * (w + GAP)) / W
            for i in range(3):
                y0 = H - T_ - (i + 1) * RH - i * RG
                axs[i, j] = fig.add_subplot(position=[x0, y0 / H, w / W, RH / H],
                                            sharex=axs[0, 0] if (i or j) else None, sharey=axs[i, 0] if j else None)
            axs[3, j] = fig.add_subplot(position=[x0, ZB / H, w / W, ZH / H], sharey=axs[3, 0] if j else None)
        x = np.arange(len(stages))
        hi = [v for q in ("zea", "human") for m in ORDER for s, _ in stages for v in trio_vals(R["se"][q][m][sv]["stages"][s])]
        mir_top = max(50, 10 * np.ceil((np.nanmax(hi) + 2) / 10))        # room for the per-trio dots
        for j, sp in enumerate(("zea", "human")):
            ac, an, am, az = axs[:, j]
            for k, m in enumerate(ORDER):
                S, off, c = R[sp]["stages"][m][sv], (k - 2) * 0.09, F.COLORS[m]
                P = R["se"][sp][m][sv]["stages"]
                kw = dict(color=c, lw=1.2, ls=F.linestyle(m), marker="o", ms=3.4, mec="black", mew=0.4, zorder=3)
                cov = [np.median(S[s]["cover"]) for s, _ in stages]
                nc = np.array([S[s]["n_child"] for s, _ in stages], float)
                mir = [nz(S[s]["mir"]) for s, _ in stages]
                ac.plot(x + off, cov, label=m, **kw)
                an.plot(x + off, nc, **kw)
                am.plot(x + off, mir, **kw)
                # the values behind each pooled value, behind the markers: each parent (cover) and each trio
                # (MIR); the call counts are totals
                for i, (s, _) in enumerate(stages):
                    trio_dots(ac, i + off, S[s]["cover"], c, s=3.5)
                    trio_dots(am, i + off, trio_vals(P[s]), c, s=3.5)
            for a in (ac, an, am):
                a.axvspan(1.6, 2.4, color="#e6e6e6", lw=0, zorder=0)
            for a in (ac, an):
                a.tick_params(axis="x", labelbottom=False)
            ac.set(ylim=(0, 100), yticks=[0, 25, 50, 75, 100])
            top = max(R[q]["stages"][m][sv][s]["n_child"] for q in ("zea", "human") for m in ORDER for s, _ in stages)
            an.set(yscale="log", ylim=(10, max(1e5, 10 ** np.ceil(np.log10(top)))))
            an.tick_params(axis="y", which="minor", left=True, right=True)
            an.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))   # no 6 pt exponents
            am.set_ylim(0, mir_top)
            am.set_xticks(x, [s for _, s in stages])
            am.set_xlim(-0.5, len(stages) - 0.5)
            ac.set_title((F.ZEA + " " if sp == "zea" else "Human ") + SVL[sv], loc="left")
            # row 4: MIR by child call length, main analysis
            labs = [b[2] for b in M.BINS]
            xs = np.arange(len(labs))
            for k, m in enumerate(ORDER):
                pts = [R[sp]["headline"][m][sv]["bins"][b] for b in labs]
                n = np.array([sum(p[1] for p in pt) for pt in pts])
                y = np.array([nz(pooled(pt)) for pt in pts])
                ok = n > 0                                            # bins without calls are left out
                xo = xs + (k - 2) * 0.09                              # callers side by side, as in rows 1-3
                az.plot(xo[ok], y[ok], color=F.COLORS[m], lw=1.2, ls=F.linestyle(m), zorder=3)
                az.scatter(xo[ok], y[ok], s=msize(n[ok]), color=tint(F.COLORS[m], 0.55), edgecolor=F.COLORS[m],
                           lw=0.6, zorder=4)
                # the circles are large, so here the per-trio dots are drawn over them
                for xx, pt, o in zip(xo, pts, ok):
                    if o:
                        trio_dots(az, xx, trio_vals(pt), F.COLORS[m], s=3.5, z=5)
            az.set_xticks(xs, ["0.1–2 kb", "2–10 kb", "10–50 kb", "≥50 kb"])
            az.set_xlim(-0.4, len(labs) - 0.6)
            az.set_xlabel("Child call length")
            if j:
                for a in axs[:, j]:
                    a.tick_params(axis="y", labelleft=False)
            for a, lab in zip(axs[:, j], "aceg" if j == 0 else "bdfh"):
                F.label(a, lab, x=-(L_ - 0.3) / w if j == 0 else -(GAP - 1.0) / w, y=1.0, va="top")
        axs[0, 0].set_ylabel("Parent genome\ncoverage (%)")
        axs[1, 0].set_ylabel("# Child calls")
        axs[2, 0].set_ylabel("MIR (%)")
        axs[3, 0].set_ylabel("MIR (%)")
        for a in axs[:, 0]:
            a.yaxis.set_label_coords(-12.6 / w, 0.5)            # right edge of the label, left of the ticks
        sz = [axs[3, 1].scatter([], [], s=msize([n]), color=F.GREY, edgecolor="black", lw=0.4) for n in (10, 100, 1000, 5000)]
        axs[3, 1].legend(sz, ["10", "100", "1,000", "5,000"], title="Child calls", loc="upper left", ncol=2,
                         handletextpad=0.3, columnspacing=0.8, labelspacing=0.9, borderpad=0.6)
        hl = dict(zip(*axs[0, 0].get_legend_handles_labels()[::-1]))
        fig.legend([hl[m] for m in ORDER], ORDER, loc="lower center", ncol=5, handlelength=2.6,
                   bbox_to_anchor=((L_ + (W - L_ - R_) / 2) / W, 0.4 / H))
        F.save(fig, "fig_filters_size_" + sv.lower())


def pick_examples(big):
    """Parental calls > 1 Mb shown in Additional file 1, Fig. S7 b-f, chosen from the headline calls by
    role: the call spanning the Tryp77 amplification, the call that is the only parental match for the most child
    duplications, a call on the Z chromosome (one copy in females), a deletion, and the longest duplication not yet
    shown."""
    out = []
    add = lambda rs: out.append(max(rs, key=lambda r: r["end"] - r["start"])) if rs else None
    pool = lambda: [r for r in big if r not in out]
    add([r for r in pool() if r["chrom"] == L.CHROM and r["start"] < L.TRYP[0] < r["end"] and r["svtype"] == "DUP"])
    rest = [r for r in pool() if r["svtype"] == "DUP"]
    if rest:
        out.append(max(rest, key=lambda r: len(r["only_child"])))
    add([r for r in pool() if r["chrom"] == "31"])
    add([r for r in pool() if r["svtype"] == "DEL"])
    while len(out) < 5 and pool():
        add([r for r in pool() if r["svtype"] == "DUP"] or pool())
    return out


def fig_oversized(R):
    big = [r for r in R["zea"]["big_calls"] if r["svtype"] == "DUP"]
    fig = plt.figure(figsize=(170 * MM, 175 * MM))
    gs = fig.add_gridspec(3, 2, hspace=0.62, wspace=0.28, left=0.07, right=0.985, top=0.96, bottom=0.06)
    ax = fig.add_subplot(gs[0, 0])
    for m, mk in (("Delly2", "o"), ("smoove", "s")):
        b = [r for r in big if r["caller"] == m]
        ax.scatter([(r["end"] - r["start"]) / 1e6 for r in b], [r["depth_cn"] for r in b], s=6, marker=mk, lw=0,
                   color=COLORS[m], alpha=0.8, zorder=3, label=f"{m} ({len(b)} calls)")
    ax.scatter(np.full(4, (L.TRYP[1] - L.TRYP[0]) / 1e6), list(R["tryp_depth_cn"].values()), s=30, marker="*",
               color=INK, zorder=4, label="Tryp77 duplication,\nresistant parents")
    ax.set(xscale="log", xlim=(0.01, 30), ylim=(0, 10), xlabel="Call length (Mb)",
           ylabel="Parent read-depth CN in call")
    ax.axhline(2, color="#BBBBBB", lw=0.6, ls="--", zorder=1)
    ax.axvline(1, color=MUTED, lw=0.6, zorder=1)
    ax.text(1.12, 9.6, "1 Mb", fontsize=6, color=MUTED, va="top")
    ax.set_xticks([0.01, 0.1, 1, 10], ["0.01", "0.1", "1", "10"])
    ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.9), handletextpad=0.3)
    ax.set_title("Parental duplication calls > 1 Mb", loc="left")
    style.label(ax, "a", x=-0.14)
    for i, r in enumerate(pick_examples(R["zea"]["big_calls"])):
        track(fig.add_subplot(gs[(i + 1) // 2, (i + 1) % 2]), r, "bcdef"[i])
    style.save(fig, "fig_oversized_calls", max_mm=None)


def track(ax, r, letter):
    rows = [w for w in L.depth_windows(r["parent"]) if w[0] == r["chrom"]]
    span = r["end"] - r["start"]
    lo, hi = max(0, r["start"] - span // 2), min(max(w[2] for w in rows), r["end"] + span // 2)
    rows = [w for w in rows if w[2] > lo and w[1] < hi]
    mid = np.array([(w[1] + w[2]) / 2 for w in rows]) / 1e6
    cn = np.array([w[3] for w in rows])
    ax.scatter(mid, np.clip(cn, 0, 6.0), s=1.2, color="#9A9A9A", lw=0, zorder=2, rasterized=True)
    ax.plot(mid, [np.median(cn[max(0, i - 15):i + 16]) for i in range(len(cn))], color=INK, lw=0.8, zorder=3)
    col = COLORS[r["caller"]]
    ax.axvspan(r["start"] / 1e6, r["end"] / 1e6, color=tint(col, 0.10), lw=0, zorder=1)
    ax.plot([r["start"] / 1e6, r["end"] / 1e6], [6.45, 6.45], color=col, lw=3, solid_capstyle="butt", zorder=4)
    kind = "duplication" if r["svtype"] == "DUP" else "deletion"
    ax.text(r["start"] / 1e6, 6.7, f"{r['caller']} {kind} call, {span / 1e6:.1f} Mb  (mean depth CN {r['depth_cn']:.2f})",
            fontsize=6, va="bottom", color=INK)
    for a, b in r["only_child"]:
        ax.plot([a / 1e6, max(b / 1e6, a / 1e6 + (hi - lo) / 1e6 * 0.004)], [0.35, 0.35], color=INK, lw=4,
                solid_capstyle="butt", zorder=4)
    if r["only_child"]:
        ax.text(r["start"] / 1e6, 0.62, "child calls matched only by this call", fontsize=5.5, color=MUTED, va="bottom")
    if r["chrom"] == L.CHROM and r["start"] < L.TRYP[0] < r["end"]:
        ax.annotate("Tryp77", xy=(L.TRYP[0] / 1e6 - span / 1e8, 5.5), xytext=(L.TRYP[0] / 1e6 - span / 1e6 * 0.2, 4.9),
                    fontsize=6, ha="center", arrowprops=dict(arrowstyle="-", lw=0.5, color=INK))
    chrom = "Z (chr. 31)" if r["chrom"] == "31" else f"chr. {r['chrom']}"
    ax.set_title(f"{r['trio'].replace('Trio', 'Trio ')} {'resistant' if r['resistant'] else 'susceptible'} parent, "
                 f"{chrom}", loc="left")
    ax.set(xlim=(lo / 1e6, hi / 1e6), ylim=(0, 7.3), yticks=range(6), xlabel="Position (Mb)", ylabel="Read-depth CN")
    ax.axhline(2, color="#BBBBBB", lw=0.6, ls="--", zorder=1)
    style.label(ax, letter, x=-0.14)


if __name__ == "__main__":
    R = json.load(open(CACHE)) if CACHE.exists() and "--recompute" not in sys.argv else compute()
    R = load_se(R)
    write_tables(R)
    for sv in M.SVTYPES:
        for f in (fig_metric_filter, fig_size_bins, fig_mir_by_size, fig_popdesign, fig_robustness, fig_species):
            f(R, sv)
    for sv in M.SVTYPES:
        fig_filter_stages(R, sv)
        fig_filters_and_size(R, sv)
    fig_oversized(R)
