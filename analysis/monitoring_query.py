"""The 12 trio moths called as new samples against panels already held, scored against ddPCR:
the growing panels and crosstrio-9 (Additional file 1, Table S15 = monitoring_fixed_query, Fig. S14)
and the nested constant-size panels (carrier titration: ten series of 11 moths, 0-8 carriers, 600 calls per
caller; Additional file 1, Section S4.2, Fig. S13)."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgb

import cnvbench.monitoring as M
from cnvbench.locus import ddpcr
from cnvbench.style import COLORS, MM, MUTED, label, save, tint

DD = ddpcr()
MOTHS = [s for _, _, s in M.TRIO_MOTHS]
CARRIER = {s: DD[s] > M.DUP for s in MOTHS}
ROLE = {s: r for _, r, s in M.TRIO_MOTHS}
fill = lambda c, a: to_rgb(c) + (a,)

# fixed query: Table S15 and Fig. S14
BGS = [("2002", "2002"), ("2002+2012", "2002+2012"), ("2002+2012+2017", "2002+2012+2017"), ("crosstrio9", "crosstrio-9")]
FQ_CALLERS = ("Delly2 merge (-y 0)", "GATK-gCNV", "CNVcaller")     # paper caller order
FQ = {c: {b: {s: M.query_cn(b, s, c) for s in MOTHS} for b, _ in BGS} for c in FQ_CALLERS}
rows = []
for c in FQ_CALLERS:
    for b, lab in BGS:
        x, y = np.array([DD[s] for s in MOTHS]), np.array([FQ[c][b][s] for s in MOTHS])
        r = "" if len(set(y)) == 1 else f"{np.corrcoef(x, y)[0, 1]:.3f}"
        rows.append((c, lab, r, f"{np.mean(np.abs(y - x)):.3f}"))
        print(f"{c:22s} {lab:15s} r = {r or '--':6s} mean |difference| = {rows[-1][3]}")
M.write_tsv("monitoring_fixed_query", ("caller", "panel", "pearson_r", "mean_abs_diff"), rows)

ROLE_COL = {"resistant": "#B2182B", "susceptible": "#2166AC", "child": "#5A5A5A"}
MARK = {"resistant": "^", "susceptible": "v", "child": "o"}
ROLE_LABEL = {"resistant": "Resistant parent", "susceptible": "Susceptible parent", "child": "Child"}
TRIOS = sorted({t for t, _, _ in M.TRIO_MOTHS})
fig, axs = plt.subplots(3, 4, figsize=(180 * MM, 122 * MM), sharex=True, sharey=True)
fig.subplots_adjust(left=0.075, right=0.99, top=0.95, bottom=0.17, hspace=0.18, wspace=0.08)
x = np.arange(len(BGS))
for i, c in enumerate(FQ_CALLERS):
    for j, tr in enumerate(TRIOS):
        ax = axs[i, j]
        for t, role, s in M.TRIO_MOTHS:
            if t == tr:
                col = ROLE_COL[role]
                ax.axhline(DD[s], color=col, ls="--", lw=0.8, alpha=0.8, zorder=1)
                ax.plot(x, [FQ[c][b][s] for b, _ in BGS], color=col, lw=1.1, marker=MARK[role], ms=4.2,
                        markerfacecolor=fill(col, 0.45), markeredgecolor=col, zorder=3)
        ax.axvspan(2.5, 3.5, color="#F0F0F0", zorder=0, lw=0)
        ax.set_ylim(0, 9.5); ax.set_yticks([0, 2, 4, 6, 8]); ax.set_xlim(-0.4, 3.4)
        ax.grid(axis="y", zorder=0); ax.set_axisbelow(True)
        if i == 0:
            ax.set_title(tr.replace("Trio", "Trio "), fontsize=7.5)
        if j == 0:
            ax.set_ylabel(f"{c}\nCalled copy number", fontsize=7)
            label(ax, "abc"[i], x=-0.42, y=1.0)
        if i == 2:
            ax.set_xticks(x, ["2002", "+2012", "+2017", "cross-\ntrio-9"], fontsize=6.2)
hs = [plt.Line2D([], [], color=ROLE_COL[r], marker=MARK[r], ms=4.2, markerfacecolor=fill(ROLE_COL[r], 0.45), lw=1.1,
                 label=ROLE_LABEL[r]) for r in MARK]
hs += [plt.Line2D([], [], color="#5A5A5A", ls="--", lw=0.8, label="ddPCR copy number"),
       plt.Rectangle((0, 0), 1, 1, color="#F0F0F0", label="Trio benchmark background")]
fig.legend(handles=hs, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.53, 0.0), fontsize=6.4,
           columnspacing=1.4, handlelength=2.0)
fig.text(0.53, 0.088, "Background panel", ha="center", fontsize=7)
save(fig, "fixed_query_stacked_y0")

# carrier titration on the nested panels: Fig. S13 and the numbers of Section S4.2
K = (0, 2, 4, 6, 8)
SERIES = range(1, 11)
TIT = {c: [(k, r, s, M.query_cn(f"nest_s{r}_k{k}", s, c)) for k in K for r in SERIES for s in MOTHS]
       for c in ("GATK-gCNV", "CNVcaller")}
OFF = {"GATK-gCNV": -0.22, "CNVcaller": 0.22}
called = lambda v: v >= M.DUP
false = lambda v: v >= M.DUP or v <= M.DEL
fig, axes = plt.subplots(1, 2, figsize=(180 * MM, 72 * MM), gridspec_kw=dict(width_ratios=[1, 1.15]))
fig.subplots_adjust(left=0.07, right=0.99, top=0.91, bottom=0.33, wspace=0.25)
xt = [f"{k}/11\n({100 * k / 11:.0f}%)" for k in K]
ax = axes[0]
for c, rows in TIT.items():
    for car, hit, ls, mk in ((True, called, "-", "o"), (False, false, "--", "s")):
        pooled, series = [], []
        for k in K:
            sel = [(r, v) for kk, r, s, v in rows if kk == k and CARRIER[s] == car]
            pooled.append(100 * np.mean([hit(v) for _, v in sel]))
            series += [(k, 100 * np.mean([hit(v) for rr, v in sel if rr == r])) for r in SERIES]
        ax.plot(K, pooled, color=COLORS[c], ls=ls, lw=1.4, marker=mk, ms=4, mfc="white" if not car else COLORS[c],
                mec=COLORS[c], zorder=3)
        ax.scatter(np.array([k for k, _ in series]) + OFF[c] * 0.5, [v for _, v in series], s=6,
                   color=tint(COLORS[c], 0.45), lw=0, zorder=2)
        print(f"{c:10s} {'carriers called' if car else 'non-carriers false'}: " + ", ".join(f"{k}: {p:.0f}%" for k, p in zip(K, pooled)))
ax.set_xticks(K, xt); ax.set_xlim(-0.6, 8.6); ax.set_ylim(-4, 104); ax.set_yticks([0, 25, 50, 75, 100])
ax.set_ylabel("Calls (%)"); ax.set_xlabel("Carriers in the 11-moth background")
ax.grid(axis="y", zorder=0); ax.set_axisbelow(True)
ax.set_title("Call outcome", loc="left"); label(ax, "a", x=-0.17)
ax = axes[1]
rng = np.random.default_rng(5)
for c, rows in TIT.items():
    for k in K:
        v = np.array([cn - DD[s] for kk, _, s, cn in rows if kk == k and CARRIER[s]])
        ax.scatter(k + OFF[c] + rng.uniform(-0.12, 0.12, len(v)), v, s=9, facecolor=fill(COLORS[c], 0.35),
                   edgecolor=COLORS[c], lw=0.5, zorder=3)
        ax.plot([k + OFF[c] - 0.17, k + OFF[c] + 0.17], [np.median(v)] * 2, color=COLORS[c], lw=2.2, zorder=4,
                solid_capstyle="butt")
        print(f"{c:10s} {k} carriers: median called - ddPCR {np.median(v):+.2f}")
    non1 = [v for kk, _, s, v in rows if kk == 8 and not CARRIER[s]]
    print(f"{c:10s} 8 carriers: {100 * np.mean([v == 1 for v in non1]):.0f}% of non-carrier calls at one copy")
ax.axhline(0, color=MUTED, lw=0.8, ls=":")
ax.set_xticks(K, xt); ax.set_xlim(-0.8, 8.8); ax.set_ylim(-7.2, 2.2); ax.set_yticks([-6, -4, -2, 0, 2])
ax.set_ylabel("Called − ddPCR copy number"); ax.set_xlabel("Carriers in the 11-moth background")
ax.grid(axis="y", zorder=0); ax.set_axisbelow(True)
ax.set_title("Carrier copy number", loc="left"); label(ax, "b", x=-0.15)
hs = [plt.Line2D([], [], color=COLORS[c], lw=1.4, marker="o", ms=4, label=c) for c in TIT]
hs += [plt.Line2D([], [], color="#555555", lw=1.4, marker="o", ms=4, label="Carriers called (CN $\\geq$ 2.5)"),
       plt.Line2D([], [], color="#555555", lw=1.4, ls="--", marker="s", ms=4, mfc="white",
                  label="Non-carriers miscalled (CN $\\geq$ 2.5 or $\\leq$ 1.5)"),
       plt.Line2D([], [], ls="", marker="o", ms=3, color="#BBBBBB", label="Single series (a) or call (b)"),
       plt.Line2D([], [], color="#555555", lw=2.2, label="Median (b)")]
fig.legend(handles=hs, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0), fontsize=6.4,
           handlelength=2.2, columnspacing=1.6)
save(fig, "gatk_carrier_titration_nested")
