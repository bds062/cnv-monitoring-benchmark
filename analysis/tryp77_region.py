"""Duplication calls around Tryp77 in the 12 H. zea trio moths, and what the filters do to them (Results,
fourth observation on ddPCR; Additional file 1, Section S3.3). Also the megabase parental calls of S3.3.

Region: the trypsin cluster on chromosome 7, Tryp079 to Tryp077 (locus.TRYP); a call counts if it overlaps it. Calls are the headline calls of each
caller (crosstrio-9; CNVkit crosstrio-3). Filters as in the main analysis (repeat filter, six-of-eight-parent
blacklist) and the earlier rule (calls recurring in all four trios). Table tryp77_region.
    python analysis/tryp77_region.py"""
from collections import defaultdict

from cnvbench import calls as C, locus as L, mir as M
from cnvbench.paths import REPEATS, TABLES
from cnvbench.style import ORDER

W0, W1 = L.TRYP
BIG = 1_000_000


def in_window(c):
    return c.chrom == L.CHROM and c.svtype == "DUP" and c.end >= W0 and c.start <= W1


def count(calls):
    return sum(in_window(c) for cs in calls.values() for c in cs)


def kb(x):
    return f"{x / 1e3:.1f}"


def gt_records(caller, background, sample):
    """(start, end, SVTYPE, GT) of the sample's non-reference chromosome 7 records (Delly2: PASS only)."""
    out = []
    for cols, r in C.vcf_records(C.vcf_path("zea", caller, background, sample)):
        if C.norm_chrom("zea", r[0]) != L.CHROM or (caller == "Delly2" and r[6] != "PASS"):
            continue
        inf, gt = C.info(r[7]), C.fmt(r)["GT"]
        if C.nonref(gt) and inf.get("SVTYPE") in ("DEL", "DUP"):
            out.append((int(r[1]), int(inf.get("END", r[1])), inf["SVTYPE"], gt))
    return out


def main():
    trios = C.trios("zea")
    role = C.roles("zea")
    carriers = {s for s, cn in L.ddpcr().items() if cn >= 2.5}
    rep = M.load_repeats("zea", REPEATS["zea"])
    raw = C.headline("zea")
    rows = []
    print(f"Duplication calls over chr7:{W0:,}-{W1:,} (Tryp079-Tryp077), 12 trio moths")
    for m in ORDER:
        r = raw[m]
        n = count(r)
        n_rep = count(M.apply_filters(r, trios, rep, bl_rule=None))
        n_bl = count(M.apply_filters(r, trios, rep, rep_frac=None))
        n_both = count(M.apply_filters(r, trios, rep))
        n_old = count(M.apply_filters(r, trios, rep, rep_frac=None, bl_rule=("trios", 4)))
        cs = [(s, c) for s, v in r.items() for c in v if in_window(c)]
        lens = sorted(c.end - c.start + 1 for _, c in cs)
        big = [c for _, c in cs if c.end - c.start + 1 > BIG]
        with_call = {s for s, _ in cs}
        per_carrier = sorted(sum(1 for s2, _ in cs if s2 == s) for s in carriers)
        rows.append((m, n, n - n_rep, n - n_bl, n_both, n - n_old, len(big),
                     kb(lens[0]) if lens else "", kb(lens[-1]) if lens else "",
                     len(with_call & carriers), len(with_call - carriers)))
        print(f"  {m:10s} {n:3d} calls ({kb(lens[0]) if lens else '-'}-{kb(lens[-1]) if lens else '-'} kb); "
              f"removed by repeat filter {n - n_rep}, by blacklist {n - n_bl}, by earlier all-four-trios rule {n - n_old}; "
              f"kept by both filters {n_both}; calls > 1 Mb {len(big)}"
              + (f" ({min(c.end - c.start + 1 for c in big) / 1e6:.1f}-{max(c.end - c.start + 1 for c in big) / 1e6:.1f} Mb)" if big else "")
              + f"; moths with a call: {len(with_call & carriers)} of {len(carriers)} carriers, {len(with_call - carriers)} non-carriers; calls per carrier {per_carrier}")
    with open(TABLES / "tryp77_region.tsv", "w") as fh:
        fh.write("caller\tcalls\trepeat_removed\tblacklist_removed\tkept_both_filters\tfour_trio_rule_removed\t"
                 "calls_over_1mb\tmin_kb\tmax_kb\tcarriers_with_call\tnoncarriers_with_call\n")
        fh.writelines("\t".join(map(str, x)) + "\n" for x in rows)

    # breakpoint callers' genotypes at the trypsin cluster (Additional file 1, Section S3.3)
    print("\nBreakpoint callers at the trypsin cluster (raw calls, genotypes)")
    for m in ("smoove", "Delly2"):
        span = defaultdict(list)                 # >1 Mb calls spanning Tryp079-Tryp077
        small = defaultdict(list)                # 5-8 kb duplications at 8.23-8.24 Mb
        for s in C.samples("zea"):
            for a, b, t, gt in gt_records(m, C.HEADLINE[m], s):
                if a < L.TRYP[0] and b > L.TRYP[1] and b - a > BIG:
                    span[s].append((a, b, t, gt))
                if t == "DUP" and 5_000 <= b - a <= 8_000 and a >= 8_230_000 and b <= 8_240_000:
                    small[s].append((a, b, gt))
        if span:
            lens = [b - a for v in span.values() for a, b, _, _ in v]
            print(f"  {m}: calls > 1 Mb spanning the cluster in {len(span)} of 12 moths "
                  f"({min(lens) / 1e6:.1f}-{max(lens) / 1e6:.1f} Mb)")
            by_call = defaultdict(list)
            for s, v in span.items():
                for a, b, t, gt in v:
                    by_call[(a, b, t)].append((s, gt))
            for (a, b, t), v in sorted(by_call.items(), key=lambda x: -len(x[1])):
                par = [s for s, _ in v if role[s] != "child"]
                print(f"    {t} {a:,}-{b:,} ({(b - a) / 1e6:.1f} Mb): {len(v)} moths, {len(par)} of 8 parents; "
                      + ", ".join(f"{role[s]}:{gt}" for s, gt in sorted(v, key=lambda x: role[x[0]])))
        sm = sorted({s for s in small})
        lens = [b - a for v in small.values() for a, b, _ in v]
        print(f"  {m}: 5-8 kb duplications at 8.230-8.240 Mb in {len(set(sm) & carriers)} of {len(carriers)} carriers and "
              f"{len(set(sm) - carriers)} non-carriers ({', '.join(role[s] for s in sorted(set(sm) - carriers))}); "
              f"{min(lens) / 1e3:.1f}-{max(lens) / 1e3:.1f} kb; genotypes "
              + ", ".join(sorted({gt for v in small.values() for *_, gt in v})))

    # megabase parental calls and their effect on MIR (Section S3.3)
    print("\nParental calls > 1 Mb after the main-analysis filters")
    for sp in ("zea", "human"):
        tr = C.trios(sp)
        head = M.apply_filters(C.headline(sp)["smoove"], tr, M.load_repeats(sp, REPEATS[sp]))
        for sv in M.SVTYPES:
            n = sum(c.svtype == sv and M.length(c) > BIG for t in tr for p in t.parents for c in head[p])
            print(f"  {C.SPECIES[sp]} smoove {sv}: {n} calls > 1 Mb in the eight parents")
        all_, cap = M.per_trio(head, tr, "DUP"), M.per_trio(head, tr, "DUP", parent_max=BIG)
        lost = sum(b[0] for b in cap) - sum(a[0] for a in all_)
        print(f"  {C.SPECIES[sp]} smoove DUP: {lost} of {sum(a[1] for a in all_)} child duplications lose their only match")
    for sp in ("zea", "human"):
        tr = C.trios(sp)
        head = M.apply_filters(C.headline(sp)["CNVkit"], tr, M.load_repeats(sp, REPEATS[sp]))
        print(f"  {C.SPECIES[sp]} CNVkit DUP MIR: {M.mir(head, tr, 'DUP'):.1f}% -> {M.mir(head, tr, 'DUP', parent_max=100_000):.1f}% "
              f"with parental calls > 100 kb kept from matching; {M.mir(head, tr, 'DUP', parent_max=BIG):.1f}% at 1 Mb")


if __name__ == "__main__":
    main()
