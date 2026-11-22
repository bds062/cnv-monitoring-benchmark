"""Mendelian inconsistency rate (MIR) of child calls: match rule, repeat filter, parent-only recurrence
blacklist, pooled MIR with its per-trio values, paired comparison of two callers within trios, size bins and
genome cover.

Calls are {sample: [Call]} for one caller (and background); matching never crosses chromosomes or types."""
import bisect
import math
from collections import defaultdict, namedtuple

import numpy as np

from .calls import norm_chrom, open_text

Match = namedtuple("Match", "metric thr min_ratio", defaults=(0.0,))
OC = Match("oc", 0.5)                           # headline: overlap coefficient >= 0.5
BLACKLIST = ("parents", 6)                       # calls in clusters carried by >= 6 of the 8 parents
REPEAT_FRAC = 0.5                                # calls more than 50% repeat-masked are removed
BL_RECIP = 0.9                                   # blacklist clusters: single linkage at >= 90% reciprocal overlap
BINS = [(100, 2_000, "0.1-2 kb"), (2_000, 10_000, "2-10 kb"), (10_000, 50_000, "10-50 kb"),
        (50_000, math.inf, ">=50 kb")]
SVTYPES = ("DUP", "DEL")


def length(c):
    return c.end - c.start + 1


def _by_key(calls):
    d = defaultdict(list)
    for i, c in enumerate(calls):
        d[(c.chrom, c.svtype)].append(i)
    return d


def matched(child, parents, rule=OC):
    """Per child call (input order): does a parental call of the same type match it under `rule`?"""
    hit = [False] * len(child)
    P = _by_key(parents)
    for key, idx in _by_key(child).items():
        plist = sorted((parents[j] for j in P.get(key, ())), key=lambda c: c.start)
        if not plist:
            continue
        bi, active = 0, []
        for i in sorted(idx, key=lambda i: child[i].start):
            a = child[i]
            while bi < len(plist) and plist[bi].start < a.end:
                active.append(plist[bi])
                bi += 1
            active = [b for b in active if b.end > a.start]
            la = a.end - a.start
            for b in active:
                lb = b.end - b.start
                inter = min(a.end, b.end) - max(a.start, b.start)
                if inter <= 0 or la <= 0 or lb <= 0:
                    continue
                score = inter / min(la, lb) if rule.metric == "oc" else min(inter / la, inter / lb)
                if score >= rule.thr and min(la, lb) / max(la, lb) >= rule.min_ratio:
                    hit[i] = True
                    break
    return hit


# ─── repeat filter ─────────────────────────────────────────────────────────────────────────────────
def load_repeats(species, path):
    """Masked intervals (BED) -> {chrom: (starts, ends, cumulative masked bp before each interval)}, merged."""
    raw = defaultdict(list)
    with open_text(path) as fh:
        for line in fh:
            p = line.split("\t")
            c = norm_chrom(species, p[0])
            if c is not None:
                raw[c].append((int(p[1]), int(p[2])))
    out = {}
    for c, iv in raw.items():
        iv.sort()
        merged = [list(iv[0])]
        for s, e in iv[1:]:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        st = [s for s, _ in merged]
        en = [e for _, e in merged]
        out[c] = (st, en, np.concatenate([[0], np.cumsum(np.subtract(en, st))]).tolist())
    return out


def _masked_before(rep, x):
    st, en, cum = rep
    i = bisect.bisect_right(st, x) - 1
    return 0 if i < 0 else cum[i] + min(x, en[i]) - st[i]


def repeat_fraction(c, repeats):
    rep = repeats.get(c.chrom)
    if rep is None or c.end <= c.start:
        return 0.0
    return (_masked_before(rep, c.end) - _masked_before(rep, c.start)) / (c.end - c.start)


# ─── recurrence blacklist ──────────────────────────────────────────────────────────────────────────
def clusters(calls):
    """Single-linkage clusters at >= BL_RECIP reciprocal overlap over all samples' calls.
    Returns (tagged [(sample, Call)], cluster root per tagged call)."""
    tagged = [(s, c) for s, cs in calls.items() for c in cs]
    parent = list(range(len(tagged)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    groups = defaultdict(list)
    for i, (_, c) in enumerate(tagged):
        groups[(c.chrom, c.svtype)].append(i)
    for idx in groups.values():
        idx.sort(key=lambda i: tagged[i][1].start)
        active = []
        for i in idx:
            a = tagged[i][1]
            # an earlier call b can reach BL_RECIP reciprocal overlap with a only while a.start - b.start
            # <= (1 - BL_RECIP) of b's length (small margin for float rounding)
            active = [j for j in active
                      if a.start - tagged[j][1].start <= 1.001 * (1 - BL_RECIP) * (tagged[j][1].end - tagged[j][1].start)]
            for j in active:
                b = tagged[j][1]
                inter = min(a.end, b.end) - max(a.start, b.start)
                la, lb = a.end - a.start, b.end - b.start
                if inter > 0 and la > 0 and lb > 0 and min(inter / la, inter / lb) >= BL_RECIP:
                    ra, rb = find(i), find(j)
                    if ra != rb:
                        parent[ra] = rb
            active.append(i)
    return tagged, [find(i) for i in range(len(tagged))]


def blacklisted(calls, trios, rule=BLACKLIST):
    """Set of (sample, Call) in clusters carried by >= k of the parents (rule ('parents', k)) or by
    >= k trios, children included (('trios', k), the earlier rule). Children never count towards the
    parents rule, but their calls in a flagged cluster are removed too."""
    kind, k = rule
    trio_of = {s: t.name for t in trios for s in (*t.parents, t.child)}
    parents = {s for t in trios for s in t.parents}
    tagged, root = clusters(calls)
    carriers = defaultdict(set)
    for (s, _), r in zip(tagged, root):
        if kind == "trios":
            carriers[r].add(trio_of[s])
        elif s in parents:
            carriers[r].add(s)
    bad = {r for r, v in carriers.items() if len(v) >= k}
    return {sc for sc, r in zip(tagged, root) if r in bad}


def apply_filters(calls, trios, repeats, rep_frac=REPEAT_FRAC, bl_rule=BLACKLIST):
    """Repeat filter and blacklist, both evaluated on the unfiltered calls; None disables one."""
    bl = blacklisted(calls, trios, bl_rule) if bl_rule else set()
    return {s: [c for c in cs if (s, c) not in bl and (rep_frac is None or repeat_fraction(c, repeats) <= rep_frac)]
            for s, cs in calls.items()}


# ─── MIR ───────────────────────────────────────────────────────────────────────────────────────────
def per_trio(calls, trios, svtype=None, rule=OC, lo=0, hi=math.inf, parent_max=None):
    """[(inconsistent, child calls)] per trio. Child calls in [lo, hi) bp; parental calls longer than
    parent_max are not allowed to match."""
    out = []
    for t in trios:
        child = [c for c in calls[t.child] if (svtype is None or c.svtype == svtype) and lo <= length(c) < hi]
        par = [c for p in t.parents for c in calls[p] if parent_max is None or length(c) <= parent_max]
        out.append((len(child) - sum(matched(child, par, rule)), len(child)))
    return out


def pct(k, n):
    return 100 * k / n if n else None


def paired(pa, pb):
    """Compare two callers on the same trios. pa, pb: [(inconsistent, child calls)] per trio, same trio order.
    Quasi-binomial logistic model logit(MIR) = trio + caller, so the family effect shared by both callers
    cancels; the dispersion is Pearson X2 / residual df (at least 1) and the caller effect is tested with t on the
    residual df (3 for four trios). Returns the trios in which b's MIR is above a's, the odds ratio of b against a
    and the two-sided p."""
    from scipy import stats
    rows = [(i, 0, k, n) for i, (k, n) in enumerate(pa)] + [(i, 1, k, n) for i, (k, n) in enumerate(pb)]
    rows = [r for r in rows if r[3]]
    T = len(pa)
    X = np.array([[float(r[0] == i) for i in range(T)] + [float(r[1])] for r in rows])
    k = np.array([r[2] for r in rows], float)
    n = np.array([r[3] for r in rows], float)
    beta = np.zeros(T + 1)
    for _ in range(100):                                  # IRLS
        mu = 1 / (1 + np.exp(-(X @ beta)))
        w = np.maximum(n * mu * (1 - mu), 1e-10)
        new = np.linalg.solve(X.T @ (w[:, None] * X), X.T @ (w * (X @ beta) + (k - n * mu)))
        if np.max(np.abs(new - beta)) < 1e-10:
            beta = new
            break
        beta = new
    mu = 1 / (1 + np.exp(-(X @ beta)))
    w = np.maximum(n * mu * (1 - mu), 1e-10)
    df = len(rows) - X.shape[1]
    phi = max(1.0, float(((k - n * mu) ** 2 / w).sum()) / df)
    se = math.sqrt(phi * np.linalg.inv(X.T @ (w[:, None] * X))[-1, -1])
    above = sum(kb / nb > ka / na for (ka, na), (kb, nb) in zip(pa, pb) if na and nb)
    return {"above": above, "trios": sum(1 for (_, na), (_, nb) in zip(pa, pb) if na and nb),
            "odds_ratio": math.exp(beta[-1]), "p": float(2 * stats.t.sf(abs(beta[-1]) / se, df))}


def summary(pt):
    inc, n = sum(x[0] for x in pt), sum(x[1] for x in pt)
    return {"mir": pct(inc, n), "inc": inc, "n": n, "per_trio": [pct(*x) for x in pt]}


def mir(calls, trios, svtype=None, **kw):
    pt = per_trio(calls, trios, svtype, **kw)
    return pct(sum(x[0] for x in pt), sum(x[1] for x in pt))


# ─── genome cover ──────────────────────────────────────────────────────────────────────────────────
def genome_size(species, fai):
    return sum(int(l.split("\t")[1]) for l in open(fai) if norm_chrom(species, l.split("\t")[0]) is not None)


def cover(calls, gsize):
    """Percentage of the genome inside the union of the calls."""
    by = defaultdict(list)
    for c in calls:
        by[c.chrom].append((c.start, c.end))
    tot = 0
    for iv in by.values():
        iv.sort()
        s, e = iv[0]
        for a, b in iv[1:]:
            if a > e:
                tot += e - s + 1
                s, e = a, b
            else:
                e = max(e, b)
        tot += e - s + 1
    return 100 * tot / gsize
