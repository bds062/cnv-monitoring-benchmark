"""Human child calls against long-read truth sets (Additional file 1, Tables S5 and S6).

Truth (GRCh38, autosomes): HGSVC2 freeze 4 deletions carried by the child, and HGSVC2 insertions annotated as
tandem duplications (TANDEM_DUP, the duplicated reference interval) as duplication truth, for the three children
with assemblies; HGSVC3 deletions (NA19240, HG02492) and Platinum Pedigree v1.2 deletions within its confident
regions (NA12878) as sensitivity checks. Calls: the headline call sets, raw and filtered, child calls >= 1 kb.
A call is true if it matches a truth variant >= 500 bp with OC >= 0.5; recall is the share of truth variants
>= 1 kb matched by a call >= 500 bp.
    python analysis/human_truthset.py"""
import bisect
from collections import defaultdict

from cnvbench import calls as C, mir as M
from cnvbench.paths import REPEATS, TABLES, TRUTH
from cnvbench.style import ORDER

SP = "human"
CHILDREN = ["NA19240", "HG02492", "NA12878"]       # the children with HGSVC2 assemblies (HG00552 has none)
MIN_CALL, MIN_MATCH, MIN_RECALL = 1_000, 500, 1_000
AUTO = {str(i) for i in range(1, 23)}


def carried(gt):
    return "1" in gt.replace("/", "|").split("|")


def read_sym_vcf(path, samples):
    """{sample: {'DEL': [Call], 'INS': {id: Call}}} of the variants each sample carries."""
    out = {s: {"DEL": [], "INS": {}} for s in samples}
    for cols, r in C.vcf_records(path):
        c = C.norm_chrom(SP, r[0])
        if c not in AUTO:
            continue
        inf = C.info(r[7])
        sv, pos, ln = inf.get("SVTYPE"), int(r[1]), abs(int(inf.get("SVLEN", 0)))
        for s in samples:
            if carried(r[cols.index(s)].split(":")[0]):
                if sv == "DEL":
                    out[s]["DEL"].append(C.Call(c, pos, pos + ln, "DEL"))
                elif sv == "INS":
                    out[s]["INS"][r[2]] = C.Call(c, pos, pos + ln, "INS")
    return out


def tandem_regions(path):
    """{insertion id: duplicated reference interval} for insertions >= 1 kb annotated TANDEM_DUP."""
    td = {}
    with C.open_text(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        iid, ilen, itd = hdr.index("ID"), hdr.index("SVLEN"), hdr.index("TANDEM_DUP")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[itd] in ("", ".", "NA") or abs(int(p[ilen])) < MIN_RECALL:
                continue
            ch, rng = p[itd].rsplit(":", 1)
            a, b = rng.split("-")
            if C.norm_chrom(SP, ch) in AUTO:
                td[p[iid]] = C.Call(C.norm_chrom(SP, ch), int(a), int(b), "DUP")
    return td


def platinum(vcf, bed):
    """Platinum Pedigree deletions of NA12878 inside its confident regions, and a containment test."""
    conf = defaultdict(list)
    with C.open_text(bed) as fh:
        for line in fh:
            p = line.split("\t")
            if C.norm_chrom(SP, p[0]) in AUTO:
                conf[C.norm_chrom(SP, p[0])].append((int(p[1]), int(p[2])))
    for v in conf.values():
        v.sort()
    starts = {c: [a for a, _ in v] for c, v in conf.items()}

    def inside(c, s, e):
        i = bisect.bisect_right(starts.get(c, []), s) - 1
        return i >= 0 and conf[c][i][1] >= e

    dels = []
    for _, r in C.vcf_records(vcf):
        c = C.norm_chrom(SP, r[0])
        if c not in AUTO or not carried(r[9].split(":")[0]):
            continue
        inf = C.info(r[7])
        if inf.get("SVTYPE") == "DEL":
            ln = abs(int(inf.get("SVLEN", 0)))
        elif len(r[3]) - len(r[4]) >= 50 and not r[4].startswith("<"):
            ln = len(r[3]) - len(r[4])
        else:
            continue
        if inside(c, int(r[1]), int(r[1]) + ln):
            dels.append(C.Call(c, int(r[1]), int(r[1]) + ln, "DEL"))
    return dels, inside


def truth_sets():
    """{(child, source, svtype): [Call >= 500 bp]} and the Platinum containment test."""
    h2 = read_sym_vcf(TRUTH / "hgsvc2" / "variants_freeze4_sv_insdel_sym.vcf.gz", CHILDREN)
    td = tandem_regions(TRUTH / "hgsvc2" / "variants_freeze4_sv_insdel.tsv.gz")
    h3 = read_sym_vcf(TRUTH / "hgsvc3" / "variants_GRCh38_sv_insdel_sym_HGSVC2024v1.0.vcf.gz", CHILDREN[:2])
    pdel, inside = platinum(TRUTH / "platinum" / "NA12878_hq_v1.2.svs.vcf.gz", TRUTH / "platinum" / "NA12878_hq_v1.2.svs.bed.gz")
    truth = {}
    for s in CHILDREN:
        truth[(s, "HGSVC2", "DEL")] = h2[s]["DEL"]
        truth[(s, "HGSVC2", "DUP")] = [td[i] for i in h2[s]["INS"] if i in td]
    for s in CHILDREN[:2]:
        truth[(s, "HGSVC3", "DEL")] = h3[s]["DEL"]
    truth[("NA12878", "Platinum", "DEL")] = pdel
    return {k: [c for c in v if c.end - c.start >= MIN_MATCH] for k, v in truth.items()}, inside


def main():
    truth, inside = truth_sets()
    trios = {t.child: t for t in C.trios(SP)}
    raw = C.headline(SP)
    rep = M.load_repeats(SP, REPEATS[SP])
    sets = {"raw": raw, "filtered": {m: M.apply_filters(raw[m], list(trios.values()), rep) for m in ORDER}}
    rows = []
    for (src, sv) in (("HGSVC2", "DEL"), ("HGSVC2", "DUP"), ("HGSVC3", "DEL"), ("Platinum", "DEL")):
        for cond, calls in sets.items():
            for m in ORDER:
                n = sup = tn = found = inc = sup_inc = sup_cons = 0
                for s in CHILDREN:
                    if (s, src, sv) not in truth:
                        continue
                    tv = truth[(s, src, sv)]
                    ch = [c for c in calls[m][s] if c.svtype == sv and c.chrom in AUTO]
                    if src == "Platinum":
                        ch = [c for c in ch if inside(c.chrom, c.start, c.end)]
                    big = [c for c in ch if M.length(c) >= MIN_CALL]
                    true = M.matched(big, [c._replace(svtype=sv) for c in tv])
                    cons = M.matched(big, [c for p in trios[s].parents for c in calls[m][p]])
                    trec = [c._replace(svtype=sv) for c in tv if c.end - c.start >= MIN_RECALL]
                    found += sum(M.matched(trec, [c for c in ch if M.length(c) >= MIN_MATCH]))
                    n, tn, sup = n + len(big), tn + len(trec), sup + sum(true)
                    inc += sum(not x for x in cons)
                    sup_inc += sum(t and not k for t, k in zip(true, cons))
                    sup_cons += sum(t and k for t, k in zip(true, cons))
                p, r = M.pct(sup, n), M.pct(found, tn)
                f1 = 2 * p * r / (p + r) / 100 if p and r else 0.0
                rows.append([src, cond, sv, m, n, p, r, f1, M.pct(inc, n), M.pct(sup_inc, inc), M.pct(sup_cons, n - inc)])
    TABLES.mkdir(parents=True, exist_ok=True)
    with open(TABLES / "truthset_human.tsv", "w") as fh:        # Tables S5 (HGSVC2) and S6
        fh.write("truth\tcalls\tsvtype\tcaller\tn\tprecision\trecall\tf1\tmir\tprecision_inconsistent\t"
                 "precision_consistent\n")
        for r in rows:
            fh.write("\t".join(map(str, r[:5])) + "\t" + "\t".join("" if x is None else f"{x:.3f}" if i == 2 else f"{x:.1f}"
                                                                 for i, x in enumerate(r[5:])) + "\n")
    print("wrote", TABLES / "truthset_human.tsv")


if __name__ == "__main__":
    main()
