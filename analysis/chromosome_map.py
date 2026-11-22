"""Chromosome correspondence between the earlier GA-R assembly (ASM2234304v1, GCA_022343045.1) and the reference
(ilHelZeax1.1, GCF_022581195.2); Additional file 1, Table S9 and Methods (chromosome 9 of the earlier assembly is
chromosome 7 of the reference). From the wfmash alignments of each earlier chromosome to the whole reference
(data/wfmash/alignments.tsv.gz): for each query-target pair, sum the aligned query and target bases; the target with
the highest query coverage is the match. Table chromosome_map.
    python analysis/chromosome_map.py"""
import collections
import csv
import gzip

from cnvbench.paths import REPO, TABLES

W = REPO / "data" / "wfmash"


def main():
    names = {r["accession"]: r["chromosome"] for r in csv.DictReader(open(W / "chromosome_names.tsv"), delimiter="\t")}
    qbp, tbp, qlen = collections.Counter(), collections.Counter(), {}
    with gzip.open(W / "alignments.tsv.gz", "rt") as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
            k = (r["query"], r["target"])
            qbp[k] += int(r["query_end"]) - int(r["query_start"])
            tbp[k] += int(r["target_end"]) - int(r["target_start"])
            qlen[r["query"]] = int(r["query_length"])
    rows = []
    for q in sorted(qlen, key=lambda a: int(names[a])):
        t = max((k for k in qbp if k[0] == q), key=lambda k: qbp[k])
        second = sorted((qbp[k] for k in qbp if k[0] == q), reverse=True)[1:2] or [0]
        rows.append(dict(query_chrom=names[q], target_chrom=names.get(t[1], t[1]), query=q, target=t[1],
                         query_length=qlen[q], mapped_query_bp=qbp[t], mapped_target_bp=tbp[t],
                         query_coverage_pct=f"{100 * qbp[t] / qlen[q]:.2f}",
                         next_best_coverage_pct=f"{100 * second[0] / qlen[q]:.2f}"))
    TABLES.mkdir(parents=True, exist_ok=True)
    with open(TABLES / "chromosome_map.tsv", "w") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"{r['query_chrom']:>3s} -> {r['target_chrom']:>3s}  {r['query_length']:>9d} {r['mapped_query_bp']:>9d} "
              f"{r['mapped_target_bp']:>9d} {r['query_coverage_pct']:>6s}%  (next best {r['next_best_coverage_pct']}%)")
    print("wrote", TABLES / "chromosome_map.tsv")


if __name__ == "__main__":
    main()
