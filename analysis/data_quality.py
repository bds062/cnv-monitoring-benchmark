"""Sequencing and alignment quality of the 24 trio genomes (Additional file 1, Table S16).

Read-level values come from `samtools stats` on 100 random 100 kb autosomal windows (seed 1; the same windows
for every sample of a species), with a duplicate estimate applied identically to both species: read pairs whose
first read shares chromosome, position, strand, mate position and insert size with another pair. Read pairs and
trimming come from the Trimmomatic logs (H. zea only; the human reads were released aligned), mapped reads and
mean aligned depth from `samtools idxstats` (reads mapped to the chromosomes x mean read length / chromosome
length).
    python analysis/data_quality.py        (samtools on PATH)"""
import random
import re
import subprocess
import tempfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from cnvbench import calls as C
from cnvbench.paths import TABLES, bam, trim_log

AUTOSOMES = {"zea": {str(i) for i in range(1, 31)}, "human": {str(i) for i in range(1, 23)}}


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout


def idxstats(path):
    return [(f[0], int(f[1]), int(f[2]), int(f[3])) for f in (l.split("\t") for l in sh(f"samtools idxstats {path}").splitlines())]


def windows(species):
    """100 random 100 kb windows on autosomes longer than 200 kb, weighted by length."""
    c = [(n, l) for n, l, _, _ in idxstats(bam(species, C.samples(species)[0]))
         if C.norm_chrom(species, n) in AUTOSOMES[species] and l > 200_000]
    random.seed(1)
    picks = random.choices(c, weights=[l for _, l in c], k=100)
    return sorted((n, a, a + 100_000) for n, a in ((n, random.randrange(l - 100_000)) for n, l in picks))


def sample_stats(species, s, bed):
    path = bam(species, s)
    sn, sec = {}, {}
    for line in sh(f"samtools view -u -M -L {bed} {path} | samtools stats -").splitlines():
        f = line.split("\t")
        if f[0] == "SN":
            sn[f[1].rstrip(":")] = float(f[2])
        elif f[0] in ("FFQ", "LFQ", "MAPQ"):
            sec.setdefault(f[0], []).append([float(x) for x in f[1:]])
    q = np.zeros(100)
    for k in ("FFQ", "LFQ"):
        for row in sec[k]:
            q[:len(row) - 1] += row[1:]                     # column i + 1 = bases of quality i
    mapq = sum(n for v, n in sec["MAPQ"] if v >= 20)
    keys = Counter()
    for line in sh(f"samtools view -M -L {bed} -f 66 -F 2304 {path}").splitlines():
        f = line.split("\t", 9)
        keys[(f[2], f[3], int(f[1]) & 48, f[7], f[8])] += 1
    idx = idxstats(path)
    chrom = [x for x in idx if C.norm_chrom(species, x[0]) is not None]
    mapped, unmapped = sum(x[2] for x in idx), sum(x[3] for x in idx)
    r = {"depth": sum(x[2] for x in chrom) * sn["average length"] / sum(x[1] for x in chrom),
         "mapped": 100 * mapped / (mapped + unmapped),
         "q30": 100 * q[30:].sum() / q.sum(), "mismatch": 100 * sn["error rate"],
         "insert": sn["insert size average"], "insert_sd": sn["insert size standard deviation"],
         "proper": sn["percentage of properly paired reads (%)"], "mapq20": 100 * mapq / sn["reads mapped"],
         "dup": 100 * (1 - len(keys) / sum(keys.values()))}
    if species == "zea":
        m = re.search(r"Input Read Pairs: (\d+) Both Surviving: (\d+)", open(trim_log(species, s)).read())
        r["pairs_M"], r["kept"] = int(m.group(1)) / 1e6, 100 * int(m.group(2)) / int(m.group(1))
    else:
        r["pairs_M"] = r["kept"] = r["mapped"] = np.nan        # released aligned, chr1-22, X and Y only
    return r


def main():
    cols = ["depth", "pairs_M", "kept", "mapped", "q30", "mismatch", "insert", "insert_sd", "proper", "mapq20", "dup"]
    TABLES.mkdir(parents=True, exist_ok=True)
    with open(TABLES / "data_quality.tsv", "w") as fh:
        fh.write("species\ttrio\trole\tsample\t" + "\t".join(cols) + "\n")
        for sp in ("zea", "human"):
            role = C.roles(sp)
            with tempfile.NamedTemporaryFile("w", suffix=".bed") as bed:
                bed.writelines(f"{n}\t{a}\t{b}\n" for n, a, b in windows(sp))
                bed.flush()
                with ThreadPoolExecutor(8) as ex:
                    res = dict(zip(C.samples(sp), ex.map(lambda s: sample_stats(sp, s, bed.name), C.samples(sp))))
            for t in C.trios(sp):
                for s in (*t.parents, t.child):
                    fh.write(f"{C.SPECIES[sp]}\t{t.name}\t{role[s]}\t{s}\t" + "\t".join(f"{res[s][c]:.2f}" for c in cols) + "\n")
    print("wrote", TABLES / "data_quality.tsv")


if __name__ == "__main__":
    main()
