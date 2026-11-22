"""Sensitivity of filtered MIR to the reciprocal overlap at which calls are clustered for the blacklist
(main analysis 0.9; Methods, Blacklist): 0.9, 0.8, 0.7 and 0.5, everything else as in the main analysis.
Table blacklist_overlap."""
from cnvbench import calls as C, mir as M
from cnvbench.paths import REPEATS, TABLES
from cnvbench.style import ORDER

rows = []
for sp in ("zea", "human"):
    trios = C.trios(sp)
    rep = M.load_repeats(sp, REPEATS[sp])
    raw = C.headline(sp)
    for m in ORDER:
        for r in (0.9, 0.8, 0.7, 0.5):
            M.BL_RECIP = r                         # read by mir.clusters() at call time
            head = M.apply_filters(raw[m], trios, rep)
            for sv in M.SVTYPES:
                n = sum(c.svtype == sv for t in trios for c in head[t.child])
                rows.append((sp, m, r, sv, n, f"{M.mir(head, trios, sv):.1f}"))
                print(*rows[-1], sep="\t", flush=True)
M.BL_RECIP = 0.9

TABLES.mkdir(parents=True, exist_ok=True)
with open(TABLES / "blacklist_overlap.tsv", "w") as fh:
    fh.write("species\tcaller\treciprocal_overlap\tsvtype\tchild_calls\tmir\n")
    fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
print("wrote", TABLES / "blacklist_overlap.tsv")
