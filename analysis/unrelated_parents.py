"""Unrelated-parent control for MIR (Additional file 1, Section S14, Table S23). Each child is scored against
its own parents (the main analysis) and against the parents of each of the other three trios of its species,
after the main-analysis filters (OC >= 0.5); calls are pooled over the 12 pairings. Also with parental calls
longer than 1 Mb kept from matching. Table unrelated_parents."""
from cnvbench import calls as C, mir as M
from cnvbench.paths import REPEATS, TABLES
from cnvbench.style import ORDER

rows = []
for sp in ("zea", "human"):
    trios = C.trios(sp)
    rep = M.load_repeats(sp, REPEATS[sp])
    raw = C.headline(sp)
    unrelated = [C.Trio(f"{t.name}x{u.name}", u.parents, t.child) for t in trios for u in trios if u is not t]
    for m in ORDER:
        head = M.apply_filters(raw[m], trios, rep)
        for sv in M.SVTYPES:
            for cap, parent_max in (("all", None), ("<=1Mb", 1_000_000)):
                own = M.mir(head, trios, sv, parent_max=parent_max)
                unrel = M.mir(head, unrelated, sv, parent_max=parent_max)
                rows.append((sp, m, sv, cap, f"{own:.1f}", f"{unrel:.1f}"))
                print(*rows[-1], sep="\t", flush=True)

TABLES.mkdir(parents=True, exist_ok=True)
with open(TABLES / "unrelated_parents.tsv", "w") as fh:
    fh.write("species\tcaller\tsvtype\tparent_calls\tmir_own_parents\tmir_unrelated_parents\n")
    fh.writelines("\t".join(r) + "\n" for r in rows)
print("wrote", TABLES / "unrelated_parents.tsv")
