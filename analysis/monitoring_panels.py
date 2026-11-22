"""The callers run on panels of the published Louisiana moths: the growing monitoring panel (Additional file 1, Fig. S26),
the other year mixes (Fig. S17) and the nested constant-size panels (Fig. 4: ten
series of 11 moths in which two non-carriers are replaced by two carriers per step, pooled over the series
at each level), each caller's copy number at the Tryp77 amplicon scored against the caller-free read-depth
truth. Tables: monitoring_moths (Table S12), monitoring_nested_panels (Table S14), monitoring_panel_calls
(carriers called and false calls in every panel, and pooled per nested level), monitoring_nested_summary
(per level: carriers called and false calls over the ten series)."""
import re

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

import cnvbench.monitoring as M
from cnvbench.locus import ddpcr
from cnvbench import figstyle as F

DEP = M.read_depth()
CARRIER = {s: cn >= M.DUP for s, (_, cn) in DEP.items()}
OK_DARK, OK_LIGHT = "#006B4F", "#C6EADC"           # correct: carrier called / non-carrier not called
BAD_DARK, BAD_LIGHT = "#D55E00", "#EE9A6A"         # incorrect: carrier missed / non-carrier falsely called
TRUE_DARK, TRUE_LIGHT = "#4D4D4D", "#D0D0D0"       # read-depth truth: carrier / non-carrier
LOOK = {"tp": (OK_DARK, OK_DARK), "fn": (BAD_DARK, BAD_DARK), "tn": (OK_LIGHT, OK_DARK),
        "fp": (BAD_LIGHT, BAD_DARK), "del": (BAD_LIGHT, BAD_DARK)}
WEDGES = ("tp", "fn", "tn", "fp")                  # false deletions are drawn with false duplications
ROWS = [("Read depth", "Read depth"), ("Delly2 merge (-y 0)", "Delly2 merge (-y 0)"), ("GATK-gCNV", "GATK-gCNV"),
        ("CNVcaller", "CNVcaller")]                       # callers in the paper's order (breakpoint, then read depth)
SEQ = [("2002", "2002"), ("2002+2012", "+ 2012"), ("2002+2012+2017", "+ 2017")]
PERM = [("2012", "2012"), ("2017", "2017"), ("2002+2017", "2002\n+ 2017"), ("2012+2017", "2012\n+ 2017"),
        ("mixed13_a", "13 from\nall (a)"), ("mixed13_b", "13 from\nall (b)")]
NEST = sorted((p for p in M.PANELS if p.startswith("nest_")), key=lambda p: tuple(map(int, re.findall(r"\d+", p))))
SERIES = sorted({int(re.findall(r"\d+", p)[0]) for p in NEST})
LEVELS = sorted({int(re.findall(r"\d+", p)[1]) for p in NEST})
CN = {p: M.panel_cn(p) for p, _ in SEQ + PERM}
CN.update({p: M.panel_cn(p) for p in NEST})
for p in CN:
    CN[p]["Read depth"] = {s: DEP[s][1] for s in M.PANELS[p]}


def pairs(p, row):
    """(moth, copy number) of a panel; nest_kK pools the ten nested panels with K carriers (moths recur)."""
    if p.startswith("nest_k"):
        return [x for r in SERIES for x in CN[f"nest_s{r}_{p[5:]}"][row].items()]
    return list(CN[p][row].items())


def counts(pairs):
    out = dict.fromkeys(("tp", "fn", "tn", "fp", "del"), 0)
    for s, v in pairs:
        out[M.outcome(CARRIER[s], v)] += 1
    return out


def pie(ax, sizes_colors):
    parts = [(n, c) for n, c in sizes_colors if n]
    ax.pie([n for n, _ in parts], colors=[c for _, c in parts], startangle=90, counterclock=False,
           wedgeprops=dict(edgecolor="white", linewidth=0.8))
    ax.add_patch(mpatches.Circle((0, 0), 1.0, fill=False, edgecolor="black", lw=F.THIN))
    ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05); ax.set_aspect("equal"); ax.axis("off")


ROW_MM, PIE_MM, TITLE_MM = 3.3, 10.5, 4.5        # one table row under the pies, pie strip, panel title (mm)
TABLE_ROWS = {"Read depth": ("Carriers",), None: ("Called", "False dup.", "False del.")}


def panel(fig, box, row, title, cols, letter, first, rng, xlabel=None):
    """One caller: a strip of pies (one per panel) with their counts as a small table, above the moths'
    copy numbers. box = (left, bottom, width, scatter height) in mm; first: leftmost column of the figure."""
    W, H = fig.get_size_inches() / F.MM
    x0, y0, w, sh = box
    ax = fig.add_axes([x0 / W, y0 / H, w / W, sh / H])
    ytop = 5.0
    yt = y0 + sh + 1.0                                   # top of the scatter's frame + gap: the table starts here
    rows = TABLE_ROWS.get(row, TABLE_ROWS[None])
    n_rows = len(TABLE_ROWS[None])
    ypie = yt + n_rows * ROW_MM + 0.6
    step = w / len(cols)
    for k, (p, _) in enumerate(cols):
        d = pairs(p, row)
        big = len(d) > 40                    # pooled nested levels: 110 moths, smaller dots
        for x, (s, v) in zip(k + rng.uniform(-0.34 if big else -0.26, 0.34 if big else 0.26, len(d)), d):
            o = M.outcome(CARRIER[s], v)
            bad = row != "Read depth" and o in ("fn", "fp", "del")
            fc, ec = ((TRUE_DARK if CARRIER[s] else TRUE_LIGHT), TRUE_DARK) if row == "Read depth" else LOOK[o]
            ax.scatter(x, min(v, ytop - 0.1), s=(9 if bad else 6) if big else (14 if bad else 11), facecolor=fc,
                       edgecolor=ec, lw=0.4 if big else 0.5, zorder=5 if bad else 3)
        n, car = len(d), sum(CARRIER[s] for s, _ in d)
        cx = x0 + (k + 0.5) * step
        sub = fig.add_axes([(cx - PIE_MM / 2) / W, ypie / H, PIE_MM / W, PIE_MM / H])
        if row == "Read depth":
            pie(sub, [(car, TRUE_DARK), (n - car, TRUE_LIGHT)])
            vals = [f"{car}/{n}"]
        else:
            c = counts(d)
            pie(sub, [(c[w_] + (c["del"] if w_ == "fp" else 0), LOOK[w_][0]) for w_ in WEDGES])
            vals = [f"{c['tp']}/{car}" if car else "–", str(c["fp"]), str(c["del"])]
        for r, v in enumerate(vals):
            fig.text(cx / W, (yt + (n_rows - r - 0.5) * ROW_MM) / H, v, ha="center", va="center", fontsize=F.SMALL)
    for r, lab in enumerate(rows):                       # row names, left of the table
        fig.text((x0 - 1.2) / W, (yt + (n_rows - r - 0.5) * ROW_MM) / H, lab, ha="right", va="center",
                 fontsize=F.SMALL)
    top = ypie + PIE_MM
    fig.text(x0 / W, (top + 1.0) / H, title, ha="left", va="bottom", fontsize=F.FS, fontweight="bold")
    fig.text((x0 - 14.0) / W, (top + 1.0) / H, letter, ha="left", va="bottom", fontsize=F.LETTER, fontweight="bold")
    for t in (M.DUP, M.DEL):
        ax.axhline(t, color="black", ls=":", lw=F.LW, zorder=1)
    ax.set_xlim(-0.5, len(cols) - 0.5); ax.set_ylim(0.6, ytop)
    ax.set_xticks(range(len(cols)), [lab for _, lab in cols])
    ax.tick_params(axis="x", top=False)
    ax.set_yticks([1, 2, 3, 4])
    ax.tick_params(labelleft=first)
    if first:
        ax.set_ylabel("Copy number")
    if xlabel:
        ax.set_xlabel(xlabel)
    return top + 1.0 + TITLE_MM - y0                     # height of the whole panel (mm)


def legend(fig, W, H):
    sw = lambda f, e: plt.Line2D([], [], ls="", marker="o", ms=5, mfc=f, mec=e, mew=0.5)
    hs = [(sw(OK_DARK, OK_DARK), "Carrier, called"), (sw(OK_LIGHT, OK_DARK), "Non-carrier, not called"),
          (sw(BAD_DARK, BAD_DARK), "Carrier, missed"), (sw(BAD_LIGHT, BAD_DARK), "Non-carrier, false dup. or del."),
          (sw(TRUE_DARK, TRUE_DARK), "Read depth: carrier ($\\geq$ 2.5 copies)"),
          (sw(TRUE_LIGHT, TRUE_DARK), "Read depth: non-carrier"),
          (plt.Line2D([], [], color="black", ls=":", lw=F.LW), "2.5 and 1.5 copies")]
    fig.legend([h for h, _ in hs], [t for _, t in hs], loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.4 / H),
               handlelength=1.6, columnspacing=1.0)


def figure(cols, name, shape, xlabel=None):
    """Panels a-d (read depth, then the callers) on a shape = (rows, columns) grid, 170 mm wide."""
    nr, nc = shape
    W, L_, R_, GAP = F.FULL, 15.5, 1.0, 15.5             # mm; the gap holds the table's row names
    LEG, XL, SH, RGAP = 11.0, 9.5 if xlabel else 5.5, 30.0, 3.0
    w = (W - L_ - R_ - (nc - 1) * GAP) / nc
    ph = SH + 1.0 + len(TABLE_ROWS[None]) * ROW_MM + 0.6 + PIE_MM + 1.0 + TITLE_MM
    H = LEG + nr * (XL + ph) + (nr - 1) * RGAP + 0.5
    rng = np.random.default_rng(3)
    with F.paper():
        fig = plt.figure(figsize=(W * F.MM, H * F.MM))
        for i, (row, title) in enumerate(ROWS):
            r, c = divmod(i, nc)
            y0 = H - 0.5 - (r + 1) * ph - r * (XL + RGAP)
            panel(fig, (L_ + c * (w + GAP), y0, w, SH), row, title, cols, "abcd"[i], c == 0, rng,
                  xlabel if r == nr - 1 else None)
        legend(fig, W, H)
        F.save(fig, name)


if __name__ == "__main__":
    figure(SEQ, "timeseries_monitoring_sequential_y0", (1, 4))
    figure(PERM, "timeseries_monitoring_permutations_y0", (2, 2))
    figure([(f"nest_k{k}", f"{k}/11") for k in LEVELS], "timeseries_titration_nested_pooled_y0", (2, 2),
           xlabel="Carriers in the 11-moth panel")

    # carriers called and false calls in every panel; nested levels also pooled over the ten series
    rows = []
    for p in [p for p, _ in SEQ + PERM] + NEST + [f"nest_k{k}" for k in LEVELS]:
        for c in ("Read depth",) + M.CALLERS:
            d = pairs(p, c)
            n = counts(d)
            rows.append((p, c, len(d), sum(CARRIER[s] for s, _ in d), n["tp"], n["fp"], n["del"]))
    M.write_tsv("monitoring_panel_calls", ("panel", "caller", "moths", "carriers", "carriers_called", "false_dup", "false_del"), rows)
    for r in rows:
        if r[1] != "Read depth":
            print(f"{r[0]:16s} {r[1]:24s} carriers called {r[4]}/{r[3]}, false dup. {r[5]}, false del. {r[6]}")


    def dup53(moth):
        """Delly2 calls the 5.3 kb tandem duplication of the trypsin cluster (PASS, non-reference, >= 80% reciprocal overlap)."""
        a0, b0 = 8_232_775, 8_238_054
        for x in M.bcf(M.MON / "chr7" / moth / f"{moth}.bcf", "%POS\t%INFO/END\t%INFO/SVTYPE\t%FILTER[\t%GT]\n"):
            a, b, t, f, g = x.split("\t")
            if t == "DUP" and f == "PASS" and g not in ("0/0", "./.") and \
                    min(int(b), b0) - max(int(a), a0) >= 0.8 * max(int(b) - int(a), b0 - a0):
                return True
        return False


    # Table S12: the 35 published moths
    DCALL = {s: dup53(s) for s in M.PANELS["2002+2012+2017"]}
    rows = [(M.YEAR[s], s, f"{DEP[s][0]:.1f}", f"{DEP[s][1]:.2f}", "yes" if CARRIER[s] else "no", f"{M.delly_alone(s):.2f}",
             "yes" if DCALL[s] else "no") for s in M.PANELS["2002+2012+2017"]]
    M.write_tsv("monitoring_moths", ("year", "sra_run", "depth", "depth_cn", "carrier", "delly2_cn", "dup_5.3kb"), rows)
    for y in M.YEARS:
        print(f"{y}: {sum(CARRIER[s] for s in M.PANELS[y])}/{len(M.PANELS[y])} carriers")
    print(f"non-carriers <= {max(cn for s, (_, cn) in DEP.items() if not CARRIER[s]):.2f}, carriers >= "
          f"{min(cn for s, (_, cn) in DEP.items() if CARRIER[s]):.2f}; depth {min(d for d, _ in DEP.values()):.1f}-"
          f"{max(d for d, _ in DEP.values()):.1f}x (mean {np.mean([d for d, _ in DEP.values()]):.1f})")
    car, non = [s for s in DEP if CARRIER[s]], [s for s in DEP if not CARRIER[s]]
    print(f"5.3 kb DUP: {sum(DCALL[s] for s in car)}/{len(car)} carriers, {sum(DCALL[s] for s in non)}/{len(non)} non-carriers, "
          f"{sum(DCALL[s] for s in M.PANELS['2002'])}/13 moths of 2002")
    dd = ddpcr()
    tc = [s for _, _, s in M.TRIO_MOTHS if dd[s] > M.DUP]
    tn = [s for _, _, s in M.TRIO_MOTHS if dd[s] <= M.DUP]
    print(f"5.3 kb DUP in the trios: {sum(map(dup53, tc))}/{len(tc)} ddPCR carriers, {sum(map(dup53, tn))}/{len(tn)} non-carriers")

    q = [x.split("\t") for s in car for x in M.bcf(M.MON / "chr7" / s / f"{s}.cnvl.bcf", "%CHROM\t%POS\t%INFO/END\t%QUAL\t%FILTER\n")]
    q = [(float(x[3]), x[4]) for x in q if M.on_amplicon(x[0], int(x[1]), int(x[2]))]
    print(f"carriers' own Delly2 CNV calls at Tryp77: quality {min(v for v, _ in q):.0f}-{max(v for v, _ in q):.0f}, "
          f"{sum(f == 'LowQual' for _, f in q)} of {len(q)} LowQual")


    def tryp77_qual(moth):
        """Highest quality of a moth's own Delly2 CNV calls (breakpoint-refined) over the ddPCR amplicon."""
        q = [x.split("\t") for x in M.bcf(M.MON / "chr7" / moth / f"{moth}.cnvl.bcf", "%CHROM\t%POS\t%INFO/END\t%QUAL\n")]
        return max((float(x[3]) for x in q if M.on_amplicon(x[0], int(x[1]), int(x[2]))), default=float("nan"))


    tq = {s: tryp77_qual(s) for s in tc}
    print("trio carriers' own Delly2 CNV calls at Tryp77 (chr 7 realignments): "
          + ", ".join(f"{s} {dd[s]:.1f} copies: {tq[s]:.0f}" for s in sorted(tc, key=lambda s: -dd[s])))

    # per nested level: carriers called (share per series: mean, min, max) and false calls (per panel)
    rows = []
    for c in M.CALLERS:
        for k in LEVELS:
            per = [counts(CN[f"nest_s{r}_k{k}"][c].items()) for r in SERIES]
            car = [n["tp"] + n["fn"] for n in per]
            pct = [100 * n["tp"] / m for n, m in zip(per, car) if m]
            false = [n["fp"] + n["del"] for n in per]
            rows.append((c, k, len(SERIES), sum(n["tp"] for n in per), sum(car),
                         *((f"{np.mean(pct):.0f}", f"{min(pct):.0f}", f"{max(pct):.0f}") if pct else ("", "", "")),
                         sum(false), f"{np.mean(false):.1f}", min(false), max(false)))
    M.write_tsv("monitoring_nested_summary", ("caller", "carriers", "series", "carriers_called", "carriers_total",
                "called_pct_mean", "called_pct_min", "called_pct_max", "false_total", "false_mean", "false_min", "false_max"), rows)

    # Table S14: the nested series (years at 0 and 8 carriers; carriers added and non-carriers removed per step)
    years = lambda p: "/".join(str(sum(M.YEAR[s] == y for s in M.PANELS[p])) for y in M.YEARS)
    rows = []
    for r in SERIES:
        P = [M.PANELS[f"nest_s{r}_k{k}"] for k in LEVELS]
        step = lambda f: "; ".join(", ".join(f(a, b)) for a, b in zip(P, P[1:]))
        rows.append((r, years(f"nest_s{r}_k{LEVELS[0]}"), years(f"nest_s{r}_k{LEVELS[-1]}"), ", ".join(sorted(P[0])),
                     step(lambda a, b: [s for s in b if s not in a]), step(lambda a, b: [s for s in a if s not in b])))
    M.write_tsv("monitoring_nested_panels", ("series", "years_k0", "years_k8", "starting_panel", "carriers_added",
                "noncarriers_removed"), rows)
