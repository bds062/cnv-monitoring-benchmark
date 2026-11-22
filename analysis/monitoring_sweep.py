"""CNVcaller frequency sweep: 16 populations of 13 moths holding one to four resistant trio parents
(Additional file 1, Section S4.1, Fig. S12, Table S13 = table monitoring_sweep_populations).
Tryp77 counts as found if a CNV region covers the whole ddPCR amplicon."""
import re

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgb

import cnvbench.monitoring as M
from cnvbench.calls import fmt, vcf_records
from cnvbench.locus import AMPLICON, ddpcr
from cnvbench.style import COLORS, MM, MUTED, label, save

DD = ddpcr()
ROLE = {s: (t, r) for t, r, s in M.TRIO_MOTHS}
POPS = sorted((p for p in M.PANELS if p.startswith("sweep_")), key=lambda p: tuple(map(int, re.findall(r"\d+", p))))


def locus(pop):
    """{member: (CP, allele types of its genotype)} at the CNV region covering the amplicon; None if there is none."""
    for cols, r in vcf_records(M.MON / "panels" / pop / "cnvcaller" / f"{pop}.vcf.gz"):
        if M.on_chr7(r[0]) and int(r[1]) <= AMPLICON[0] and M.end(r) >= AMPLICON[1]:
            alt = ["ref"] + ["del" if a == "CN0" else "gain" for a in r[4].split(",")]
            return {n: (float(fmt(r, 9 + i)["CP"]), {alt[int(a)] for a in re.split("[/|]", fmt(r, 9 + i)["GT"]) if a != "."})
                    for i, n in enumerate(cols[9:])}
    return None


RES = {}
for p in POPS:
    k, rep = map(int, re.findall(r"\d+", p))
    car = [s for s in M.PANELS[p] if ROLE.get(s, ("", ""))[1] == "resistant"]
    RES[p] = dict(k=k, rep=rep, carriers=car, non=[s for s in M.PANELS[p] if s not in car], locus=locus(p))
    assert len(car) == k

print("carriers  found  carriers duplicated  CP/ddPCR  non-carriers non-reference")
for k in (1, 2, 3, 4):
    rs = [r for r in RES.values() if r["k"] == k]
    found = [r for r in rs if r["locus"]]
    dup = sum("gain" in r["locus"][s][1] for r in found for s in r["carriers"])
    ratio = [r["locus"][s][0] / DD[s] for r in found for s in r["carriers"]]
    nonref = sum(r["locus"][s][1] != {"ref"} for r in found for s in r["non"])
    print(f"{k}/13 ({100 * k / 13:.1f}%)  {len(found)}/{len(rs)}  {dup}/{k * len(rs)}  "
          f"{np.mean(ratio) if ratio else float('nan'):.2f}  {nonref}/{sum(len(r['non']) for r in found)}")
alone = {s: M.delly_alone(s) for r in RES.values() for s in r["carriers"]}
print("Delly2 alone, resistant parents: " + ", ".join(f"{s} {v:.2f}" for s, v in sorted(alone.items())))
M.write_tsv("monitoring_sweep_populations", ("carriers", "replicate", "resistant_parents", "n2002", "susceptible_parents"),
            [(r["k"], r["rep"], ", ".join(ROLE[s][0] for s in r["carriers"]), sum(s in M.YEAR for s in r["non"]),
              sum(s in ROLE for s in r["non"])) for r in RES.values()])

# figure
C_CNV, C_DL, NONC, GRID = COLORS["CNVcaller"], COLORS["Delly2"], "#9E9E9E", "#E3E3E3"
fig, (axA, axB) = plt.subplots(1, 2, figsize=(170 * MM, 72 * MM), gridspec_kw={"width_ratios": [1.7, 1]})
rng = np.random.default_rng(1)
G, COL = 4.0, 0.6                    # level centres at k * G; the four replicate populations around each
centres = [k * G for k in (1, 2, 3, 4)]
for k, c in zip((1, 2, 3, 4), centres):
    rs = sorted((r for r in RES.values() if r["k"] == k), key=lambda r: r["rep"])
    for i, r in enumerate(rs):
        if r["locus"]:
            xc = c + (i - (len(rs) - 1) / 2) * COL
            nc = [r["locus"][s][0] for s in r["non"]]
            ca = [r["locus"][s][0] for s in r["carriers"]]
            axA.scatter(xc + rng.uniform(-0.13, 0.13, len(nc)), nc, s=9, color=NONC, lw=0, zorder=2)
            axA.scatter(xc + rng.uniform(-0.08, 0.08, len(ca)), ca, s=22, color=C_CNV, edgecolor="white", lw=0.5, zorder=3)
    if not any(r["locus"] for r in rs):
        axA.add_patch(plt.Rectangle((c - 1.1, 0.3), 2.2, 8.4, color="#F2F2F2", lw=0, zorder=0))
        axA.text(c, 4.5, f"not\ndetected\n(0 of {len(rs)})", ha="center", va="center", fontsize=6.3, color=MUTED, linespacing=1.3)
    cs = sorted({s for r in rs for s in r["carriers"]})
    n_call = sum(alone[s] >= M.DUP for r in rs for s in r["carriers"])
    axA.scatter(np.full(len(cs), c + 1.45), [alone[s] for s in cs], marker="D", s=14, color=C_DL, edgecolor="white", lw=0.4, zorder=3)
    axA.text(c, 0.35, f"Delly2: {n_call}/{sum(len(r['carriers']) for r in rs)}", ha="center", va="bottom", fontsize=5.8, color=C_DL)
xmin, xmax = centres[0] - 1.6, centres[-1] + 1.9
axA.axhline(2, color="#BDBDBD", lw=0.7, zorder=1)
axA.text(xmin + 0.15, 2.1, "diploid", fontsize=5.8, color=MUTED, ha="left", va="bottom", zorder=3)
axA.axvline(0.1 * 13 * G, color=MUTED, lw=0.7, ls="--", zorder=1)
axA.text(0.1 * 13 * G + 0.15, 1.05, "10% (-f 0.1)", fontsize=6, color=MUTED, va="bottom")
axA.set_xticks(centres, [f"{k} of 13\n({k / 13 * 100:.1f}%)" for k in (1, 2, 3, 4)])
axA.tick_params(axis="x", length=0)
axA.set_xlim(xmin, xmax); axA.set_ylim(0, 9); axA.set_yticks([0, 2, 4, 6, 8])
axA.set_ylabel("CNVcaller copy number (CP) at Tryp077")
axA.set_xlabel("Resistant adults in a population of 13 (4 replicate populations each)")
axA.grid(axis="y", color=GRID, lw=0.4, zorder=0); axA.set_axisbelow(True)
axA.scatter([], [], marker="D", s=14, color=C_DL, edgecolor="white", lw=0.4, label="Delly2 (carrier, alone)")
axA.scatter([], [], s=22, color=C_CNV, edgecolor="white", lw=0.5, label="Resistant parent (carrier)")
axA.scatter([], [], s=9, color=NONC, lw=0, label="Non-carrier")
axA.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 1.02), ncol=3, handletextpad=0.3, columnspacing=1.0,
           fontsize=7)
axA.set_title("Genotypes at Tryp077 in each population", loc="left")
a, b = np.array([(DD[s], r["locus"][s][0]) for r in RES.values() if r["locus"] for s in r["carriers"]]).T
axB.plot([2, 9.5], [2, 9.5], ls=":", color=MUTED, lw=0.8, zorder=1)
axB.text(6.8, 7.6, "y = x", fontsize=6, color=MUTED, ha="right", va="bottom")
axB.scatter(a, b, s=24, facecolor=to_rgb(C_CNV) + (0.3,), edgecolor=to_rgb(C_CNV) + (0.9,), lw=0.9, zorder=3)
axB.set_xlim(2, 9.5); axB.set_ylim(2, 9.5); axB.set_xticks([2, 4, 6, 8]); axB.set_yticks([2, 4, 6, 8])
axB.set_aspect("equal"); axB.grid(color=GRID, lw=0.4, zorder=0); axB.set_axisbelow(True)
axB.set_xlabel("ddPCR copy number"); axB.set_ylabel("CNVcaller CP")
axB.set_title(f"Carrier copy number vs ddPCR (n = {len(a)})", loc="left")
fig.tight_layout(w_pad=2.5)
label(axA, "a", x=-0.1, y=1.02); label(axB, "b", x=-0.24, y=1.02)
save(fig, "fig_freq_detect_published")
