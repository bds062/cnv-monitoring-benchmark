"""A simple monitoring rule on chromosome 7 of the 35 published moths: every CNV site of >= 100 bp a
caller reports (gains and losses separately) is flagged if its carrier share rises with year over 2002,
2012 and 2017, flagged sites are ranked by an exact one-sided Cochran-Armitage trend test, and sites at a
Benjamini-Hochberg FDR below 5% (over all of the caller's sites) are reported; the headline is the rank of
the best Tryp77 site (Fig. 5; Delly2 run three ways: Additional file 1, Fig. S18). Table monitoring_scan;
table monitoring_scan_drift repeats the ranking with a drift-adjusted test (sensitivity check; Table S22)."""
from collections import defaultdict
from math import comb

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.legend_handler import HandlerTuple
from matplotlib.ticker import MaxNLocator
from scipy import stats

import cnvbench.monitoring as M
from cnvbench import figstyle
from cnvbench.calls import fmt, vcf_records

P = "2002+2012+2017"
D = M.MON / "panels" / P
MOTHS = M.PANELS[P]
N = {y: sum(M.YEAR[s] == y for s in MOTHS) for y in M.YEARS}
MIN_LEN = 100
GAIN_LOSS = {"gain": lambda c: c >= M.DUP, "loss": lambda c: c <= M.DEL}


def cluster(calls):
    """Per-moth calls [(start, end, kind, moth)] joined into sites: same kind, >= 50% reciprocal overlap,
    single linkage; a site spans the medians of its calls' boundaries."""
    out = []
    for kind in ("gain", "loss"):
        cs = sorted(c for c in calls if c[2] == kind)
        par = list(range(len(cs)))

        def find(i):
            while par[i] != i:
                par[i] = par[par[i]]; i = par[i]
            return i
        active = []
        for i, (a, b, _, _) in enumerate(cs):
            active = [j for j in active if cs[j][1] >= a]
            for j in active:
                a2, b2 = cs[j][:2]
                ov = min(b, b2) - max(a, a2) + 1
                if ov > 0 and ov >= 0.5 * (b - a + 1) and ov >= 0.5 * (b2 - a2 + 1):
                    par[find(i)] = find(j)
            active.append(i)
        grp = defaultdict(list)
        for i, c in enumerate(cs):
            grp[find(i)].append(c)
        out += [(int(np.median([c[0] for c in g])), int(np.median([c[1] for c in g])), kind, {c[3] for c in g})
                for g in grp.values()]
    return out


def gatk():
    calls = []
    for s in MOTHS:
        for _, r in vcf_records(D / "gatk" / f"{s}.segments.vcf.gz"):
            a, b, cn = int(r[1]), M.end(r), int(fmt(r)["CN"])
            if M.on_chr7(r[0]) and cn != 2 and b - a + 1 >= MIN_LEN:
                calls.append((a, b, "gain" if cn > 2 else "loss", s))
    return cluster(calls)


def cnvcaller():
    out = []
    for cols, r in vcf_records(D / "cnvcaller" / f"{P}.vcf.gz"):
        a, b = int(r[1]), M.end(r)
        if M.on_chr7(r[0]) and b - a + 1 >= MIN_LEN:
            cp = {n: float(fmt(r, 9 + i)["CP"]) for i, n in enumerate(cols[9:])}
            out += [(a, b, k, car) for k, hit in GAIN_LOSS.items() if (car := {n for n, x in cp.items() if hit(x)})]
    return out


def delly_per_moth():
    calls = []
    for s in MOTHS:
        for x in M.bcf(M.MON / "chr7" / s / f"{s}.bcf", "%POS\t%INFO/END\t%INFO/SVTYPE[\t%GT]\n",
                       "-i", 'FILTER="PASS" && (INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP")'):
            a, b, t, gt = x.split("\t")
            if gt not in ("0/0", "./.", ".") and int(b) - int(a) + 1 >= MIN_LEN:
                calls.append((int(a), int(b), "gain" if t == "DUP" else "loss", s))
    return cluster(calls)


def delly_merge(minqual):
    f = D / f"delly_y{minqual}" / "filtered.bcf"
    names, out = M.bcf(f, "", "-l"), []
    for x in M.bcf(f, "%POS\t%INFO/END[\t%CN]\n", "-i", 'FILTER="PASS"'):
        p = x.split("\t")
        a, b, cn = int(p[0]), int(p[1]), [int(c) if c != "." else 2 for c in p[2:]]
        if b - a + 1 >= MIN_LEN:
            out += [(a, b, k, car) for k, hit in (("gain", lambda c: c >= 3), ("loss", lambda c: c <= 1))
                    if (car := {s for s, c in zip(names, cn) if hit(c)})]
    return out


def trend_p(k):
    """Exact one-sided (increasing) Cochran-Armitage trend test over the years (scores 0, 1, 2): given the
    moths per year and the total carriers, the carriers' year counts are multivariate hypergeometric."""
    n, x = [N[y] for y in M.YEARS], [k[y] for y in M.YEARS]
    T, obs = sum(x), sum(i * v for i, v in enumerate(x))
    p = sum(comb(n[0], a) * comb(n[1], b) * comb(n[2], T - a - b)
            for a in range(min(T, n[0]) + 1) for b in range(min(T - a, n[1]) + 1)
            if T - a - b <= n[2] and b + 2 * (T - a - b) >= obs)
    return p / comb(sum(n), T)


def bh(p):
    """Benjamini-Hochberg q-values."""
    p = np.asarray(p, float); m = len(p); o = np.argsort(p); q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def trend_p_drift(k, F):
    """Asymptotic one-sided Cochran-Armitage test whose variance adds genetic drift between sampling years:
    each year's carrier share varies by F p(1 - p) around the pooled share p, the over-dispersion adjustment
    of Spitzer et al. (2020). F = 0 is the binomial (sampling-only) test. Used only as a sensitivity check."""
    n, x = (np.array([v[y] for y in M.YEARS], float) for v in (N, k))
    c = np.arange(3) - (np.arange(3) * n).sum() / n.sum()
    p = x.sum() / n.sum()
    v = (c ** 2 * (n + n ** 2 * F) * p * (1 - p)).sum()
    return 1.0 if v == 0 else float(stats.norm.sf((c * x).sum() / np.sqrt(v)))


def score(sites):
    rows = []
    for a, b, kind, car in sites:
        k = {y: sum(M.YEAR[s] == y for s in car) for y in M.YEARS}
        x, n = np.array([k[y] for y in M.YEARS]), np.array([N[y] for y in M.YEARS])
        rows.append({"start": a, "end": b, "kind": kind, "k": k, "delta": k["2017"] / N["2017"] - k["2002"] / N["2002"],
                     "rise": (np.arange(3) * (x - n * x.sum() / n.sum())).sum() > 0, "p": trend_p(k),
                     "known": kind == "gain" and min(b, M.TRYP77[1]) - max(a, M.TRYP77[0]) + 1 >= 0.5 * (b - a + 1)})
    for r, q in zip(rows, bh([r["p"] for r in rows])):
        r["q"] = q
    return rows


SITES = {"GATK-gCNV": gatk(), "CNVcaller": cnvcaller(), "Delly2 (per moth)": delly_per_moth(),
         "Delly2 merge (default)": delly_merge(200), "Delly2 merge (-y 0)": delly_merge(0)}
R, table = {}, []
for m, sites in SITES.items():
    rows = score(sites)
    flag = [r for r in rows if r["rise"]]
    sig = [r for r in flag if r["q"] < 0.05]
    order = sorted(rows, key=lambda r: (r["p"], -r["delta"]))
    rank = next((i + 1 for i, r in enumerate(order) if r["known"]), "")
    R[m] = dict(sites=rows, n=len(rows), rank=rank, flag=(len(flag), sum(r["known"] for r in flag)),
                sig=(len(sig), sum(r["known"] for r in sig)))
    table.append((m, len(rows), sum(r["known"] for r in rows), *R[m]["flag"], *R[m]["sig"], rank))
    print(f"{m:24s} sites {len(rows):5d} (Tryp77 {table[-1][2]}) | rising {len(flag):4d} (Tryp77 {R[m]['flag'][1]}) | "
          f"FDR < 5%: {len(sig):3d} (Tryp77 {R[m]['sig'][1]}) | best Tryp77 rank {rank or '-'}")
M.write_tsv("monitoring_scan", ("caller", "sites", "tryp77_sites", "flagged", "flagged_tryp77", "fdr05", "fdr05_tryp77",
                                "tryp77_best_rank"), table)

# Drift sensitivity: Taylor et al. (2021) report temporal FST < 0.004 between these collections, so drift
# adds about n F = 5% to the sampling variance of a carrier share among 11-13 moths. Rerank and re-test
# with the drift-adjusted test at F = 0 (normal approximation only), 0.004 and 0.01.
drift = []
for m in ("Delly2 merge (-y 0)", "GATK-gCNV", "CNVcaller"):
    rows = R[m]["sites"]
    for test, pf in (("exact", lambda r: r["p"]),
                     *((f"asymptotic F={F}", lambda r, F=F: trend_p_drift(r["k"], F)) for F in (0, 0.004, 0.01))):
        ps = [pf(r) for r in rows]; q = bh(ps)
        order = sorted(range(len(rows)), key=lambda i: (ps[i], -rows[i]["delta"]))
        rank = next(j + 1 for j, i in enumerate(order) if rows[i]["known"])
        sig = [i for i in range(len(rows)) if rows[i]["rise"] and q[i] < 0.05]
        drift.append((m, test, rank, len(sig), sum(rows[i]["known"] for i in sig)))
        print(f"{m:24s} {test:22s} best Tryp77 rank {rank:3d} | FDR < 5%: {len(sig):3d} (Tryp77 {drift[-1][-1]})")
M.write_tsv("monitoring_scan_drift", ("caller", "test", "tryp77_best_rank", "fdr05", "fdr05_tryp77"), drift)

KNOWN_C, NONE_C = "#006B4F", "#BDBDBD"
TRYP = figstyle.TRYP
col = lambda m: figstyle.COLORS["Delly2" if m.startswith("Delly2") else m]


def figure(tools, name):
    """Fig. 5 layout (print style, 170 mm): sites on chromosome 7 per caller (a-c), flagged sites (d, e)."""
    F = figstyle
    W, H = F.FULL, 112.0
    L_, R_, GAP = 19.0, 1.0, 5.0                                   # mm
    w = (W - L_ - R_ - 2 * GAP) / 3
    TT, TH = 6.0, 44.0                                             # top row: top margin, height
    BL, BW, BG, BB, BH = 31.0, 63.0, 6.0, 19.0, 28.0               # bottom row: left, width, gap, bottom, height
    rng = np.random.default_rng(2)
    with F.paper():
        fig = plt.figure(figsize=(W * F.MM, H * F.MM))
        for i, m in enumerate(tools):
            ax = fig.add_axes([(L_ + i * (w + GAP)) / W, (H - TT - TH) / H, w / W, TH / H])
            sites = R[m]["sites"]
            for grp, c, size, mk, z in ((lambda r: not r["rise"] and not r["known"], NONE_C, 4, "o", 2),
                                        (lambda r: r["rise"] and not r["known"], col(m), 5, "o", 3),
                                        (lambda r: r["known"], KNOWN_C, 18, "D", 5)):
                rs = [r for r in sites if grp(r)]
                if rs:
                    ax.scatter([(r["start"] + r["end"]) / 2e6 for r in rs],
                               100 * np.array([r["delta"] for r in rs]) + rng.uniform(-1, 1, len(rs)),
                               s=size, color=c, marker=mk, lw=0.3 if mk == "D" else 0, edgecolor="black", zorder=z)
            known = [r for r in sites if r["known"]]
            if known:
                r = max(known, key=lambda r: r["delta"])
                ax.annotate(TRYP, ((r["start"] + r["end"]) / 2e6, 100 * r["delta"]), xytext=(5.6, 88), ha="right",
                            va="center", fontsize=F.SMALL, color=KNOWN_C,
                            arrowprops=dict(arrowstyle="-", color=KNOWN_C, lw=0.6, shrinkB=2))
            ax.axhline(0, color="black", lw=F.THIN, zorder=1)
            ax.set_xlim(0, 13.65); ax.set_ylim(-105, 105)
            ax.set_xticks([0, 2.5, 5, 7.5, 10, 12.5], ["0", "2.5", "5", "7.5", "10", "12.5"])
            ax.set_yticks([-100, -50, 0, 50, 100], ["−100", "−50", "0", "+50", "+100"])
            ax.tick_params(labelleft=i == 0)
            if i == 0:
                ax.set_ylabel("Change in carrier\nfrequency (% points)")
            ax.set_title(m, loc="left")
            ax.text(0.97, 0.96, f"{R[m]['n']:,} sites", transform=ax.transAxes, ha="right", va="top", fontsize=F.SMALL)
            rk = R[m]["rank"]
            ax.text(0.03, 0.04, f"{TRYP} rank {rk}" if rk else f"{TRYP} not reported", transform=ax.transAxes,
                    fontsize=F.SMALL, color=KNOWN_C)
            F.label(ax, "abc"[i], x=-(L_ - 0.3) / w if i == 0 else -(GAP - 0.6) / w, y=1.0)
        fig.text((L_ + (W - L_ - R_) / 2) / W, (H - TT - TH - 8.0) / H, "Chromosome 7 position (Mb)", ha="center",
                 va="center", fontsize=F.FS)
        y = np.arange(len(tools))[::-1]
        for j, (key, xl) in enumerate((("flag", "Rising sites"), ("sig", "Rising sites, FDR < 5%"))):
            ax = fig.add_axes([(BL + j * (BW + BG)) / W, BB / H, BW / W, BH / H])
            k = np.array([R[m][key][1] for m in tools]); o = np.array([R[m][key][0] for m in tools]) - k
            ax.barh(y, k, color=KNOWN_C, height=0.6, zorder=3, **F.EDGE)
            ax.barh(y, o, left=k, color=[col(m) for m in tools], height=0.6, zorder=3, **F.EDGE)
            mx = (k + o).max()
            for yi, a, b in zip(y, k, o):
                ax.text(a + b + mx * 0.025, yi, f"{a}/{a + b:,}", va="center", fontsize=F.SMALL)
            ax.set_yticks(y, tools if j == 0 else [""] * len(tools))
            ax.set_ylim(-0.6, len(tools) - 0.4)
            ax.set_xlim(0, mx * 1.25); ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=4))
            ax.grid(False); ax.grid(True, axis="x")
            ax.set_xlabel(xl)
            F.label(ax, "de"[j], x=-(BL - 0.3) / BW if j == 0 else -(BG - 0.8) / BW, y=1.0)
        hs = [plt.Line2D([], [], ls="", marker="D", ms=4.2, mfc=KNOWN_C, mec="black", mew=0.3),
              tuple(plt.Line2D([], [], ls="", marker="o", ms=3.6, color=c) for c in dict.fromkeys(map(col, tools))),
              plt.Line2D([], [], ls="", marker="o", ms=3.6, color=NONE_C)]
        fig.legend(hs, [f"{TRYP} duplication", "Other site, rising", "Other site, not rising"],
                   handler_map={tuple: HandlerTuple(ndivide=None, pad=0.2)}, loc="lower center", ncol=3,
                   bbox_to_anchor=((L_ + (W - L_ - R_) / 2) / W, 0.4 / H), handlelength=2.2)
        F.save(fig, name)


if __name__ == "__main__":
    figure(["Delly2 merge (-y 0)", "GATK-gCNV", "CNVcaller"], "monitoring_fp_headline")   # paper caller order
    figure(["Delly2 (per moth)", "Delly2 merge (default)", "Delly2 merge (-y 0)"], "monitoring_fp_delly")
