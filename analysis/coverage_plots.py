"""Read depth, GC content and mapping quality around Tryp77 for each H. zea trio, on the earlier GA-R assembly
(GCA_022343045.1, chromosome 9, 5-6 Mb) and on the reference used here (GCF_022581195.2, chromosome 7,
7.5-8.6 Mb, drawn in the orientation of the earlier assembly): Trio1_plot ... Trio4_plot (Additional file 1).

10 kb windows: reads overlapping the window, depth relative to the mean of the region's windows outside the
amplified block, GC fraction of the reference, and mean mapping quality of the reads.
    python analysis/coverage_plots.py"""
import numpy as np
import matplotlib.pyplot as plt
import pysam

from cnvbench import calls as C, style
from cnvbench.locus import bam_contig, norm_chrom
from cnvbench.paths import REF, bam
from cnvbench.style import MM

W = 10_000
# reference, BAM kind, contig (None: chromosome 7 found by name), region, block excluded from the baseline
REFS = [("GCA_022343045.1, chr. 9", "zea_gca", "bam_gca", "CP089896.1", (5_000_000, 6_000_000), (5_200_000, 5_400_000)),
        ("GCF_022581195.2, chr. 7", "zea", "bam", None, (7_500_000, 8_600_000), (8_200_000, 8_400_000))]
ROWS = ["Read depth (reads)", "Relative read depth", "GC fraction", "Mean MAPQ"]
COLOR = {"resistant": "#1565C0", "susceptible": "#2E7D32", "child": "#C62828"}
LABEL = {"resistant": "Resistant parent", "susceptible": "Susceptible parent", "child": "Offspring"}


def tracks(path, fasta, contig, region, block):
    """(window starts, reads, relative depth, GC, mean MAPQ) over the region."""
    starts = np.arange(region[0], region[1] - W + 1, W)
    with pysam.FastaFile(str(fasta)) as fa:               # the FASTA keeps accession names; the BAMs may not
        fctg = contig if contig in fa.references else next(c for c in fa.references
                                                           if norm_chrom("zea", c) == norm_chrom("zea", contig))
        seq = fa.fetch(fctg, region[0], region[1]).upper()
    gc = np.array([(seq[a - region[0]:a - region[0] + W].count("G") + seq[a - region[0]:a - region[0] + W].count("C")) / W
                   for a in starts])
    reads, mapq = [], []
    with pysam.AlignmentFile(str(path)) as fh:
        for a in starts:
            q = [r.mapping_quality for r in fh.fetch(contig, a, a + W) if not r.is_unmapped]
            reads.append(len(q))
            mapq.append(np.mean(q) if q else 0)
    reads = np.array(reads, dtype=float)
    out = (starts < block[0]) | (starts + W > block[1])
    return starts, reads, reads / reads[out].mean(), gc, np.array(mapq)


def plot(trio):
    role = C.roles("zea")
    fig, axes = plt.subplots(4, 2, figsize=(170 * MM, 200 * MM), sharey="row")
    for j, (title, ref, kind, contig, region, block) in enumerate(REFS):
        for s in (trio.child, *trio.parents):
            path = bam("zea", s, kind)
            ctg = contig or bam_contig(path, "zea", "7")
            x, *vals = tracks(path, REF[ref], ctg, region, block)
            for ax, v in zip(axes[:, j], vals):
                ax.plot((x + W / 2) / 1e6, v, color=COLOR[role[s]], lw=0.8, alpha=0.85, label=LABEL[role[s]])
        for i, ax in enumerate(axes[:, j]):
            ax.grid(True)
            ax.set_title(f"{ROWS[i]}, {title}", loc="left", fontsize=7)
            if j == 1:
                ax.invert_xaxis()                        # GCF chromosome 7 runs opposite to GCA chromosome 9 here
        axes[-1, j].set_xlabel("Position (Mb)")
    for i, ax in enumerate(axes[:, 0]):
        ax.set_ylabel(ROWS[i].split(" (")[0])
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="lower center", ncol=3, frameon=False)
    style.save(fig, f"{trio.name}_plot")


if __name__ == "__main__":
    for t in C.trios("zea"):
        plot(t)
