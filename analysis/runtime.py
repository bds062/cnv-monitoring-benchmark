"""Run time and memory per caller step (Additional file 1, Tables S19 and S20; Results, run time and memory).
Trio benchmark: the fixed-configuration re-timing in data/runtime/fixed_cores.tsv (8 cores, one CPU type, three
repeats; workflow/runtime/timed.sh), summarised as wall-clock median (min-max), median CPU hours and the largest
peak memory -> results/tables/runtime_fixed.tsv. Preparation of the published moths and the monitoring analyses:
the as-run per-job records in data/runtime/steps.tsv, summarised as wall-clock median (min-max), cores allocated and
median / maximum peak memory -> results/tables/runtime.tsv. Where peaks were recorded per job rather than per unit,
the record's note gives their median and maximum, which are reported instead.
    python analysis/runtime.py"""
import csv
import re
import statistics as st

from cnvbench.paths import REPO, TABLES

STEPS = REPO / "data" / "runtime" / "steps.tsv"
FIXED = REPO / "data" / "runtime" / "fixed_cores.tsv"
PEAK = re.compile(r"jobs: median ([\d.]+) GB, max ([\d.]+) GB")


def load(path=STEPS):
    return list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t"))


def fmt_h(h):
    return "" if h == "" else f"{h * 60:.0f} min" if h < 1 else f"{h:.1f} h"


def summarise(recs):
    groups = {}
    for r in recs:
        groups.setdefault((r["species"], r["caller"], r["step"], r["unit"]), []).append(r)
    out = []
    for (sp, caller, step, unit), rs in groups.items():
        wall = [float(r["wall_h"]) for r in rs if r["wall_h"]]
        rss = [float(r["rss_gb"]) for r in rs if r["rss_gb"]]
        note = rs[0]["note"]
        m = PEAK.search(note)
        row = dict(species=sp, caller=caller, step=step, unit=unit, n=len(wall), source=rs[0]["source"],
                   cores=rs[0]["cores"], note=note)
        if wall:
            row.update(wall_med_h=st.median(wall), wall_min_h=min(wall), wall_max_h=max(wall))
        if rss:
            row.update(rss_med_gb=st.median(rss), rss_max_gb=max(rss))
        elif m:
            row.update(rss_med_gb=float(m.group(1)), rss_max_gb=float(m.group(2)))
        out.append(row)
    return out


def wall_s(w):
    """GNU time elapsed (h:mm:ss or m:ss.ss) in seconds."""
    p = [float(x) for x in w.split(":")]
    return p[0] * 3600 + p[1] * 60 + p[2] if len(p) == 3 else p[0] * 60 + p[1]


def summarise_fixed(path=FIXED):
    """One row per species, step and label (GATK-gCNV human shards pooled; postprocessing of the three trio
    members pooled): n, wall median/min/max (h), median CPU hours, largest peak RSS (GB)."""
    groups = {}
    for r in load(path):
        if r["rc"] != "0":
            continue
        step = "gatk_shard" if r["step"].startswith("gatk_shard") else r["step"]
        label = "post" if r["label"].startswith("post_") else r["label"]
        groups.setdefault((r["species"], step, label), []).append(r)
    out = []
    for (sp, step, label), rs in sorted(groups.items()):
        wall = [wall_s(r["wall"]) / 3600 for r in rs]
        out.append(dict(species=sp, step=step, label=label, n=len(rs), wall_med_h=st.median(wall),
                        wall_min_h=min(wall), wall_max_h=max(wall), cores=rs[0]["cores"],
                        cpu_med_h=st.median(float(r["cpu_s"]) for r in rs) / 3600,
                        rss_max_gb=max(int(r["maxrss_kb"]) for r in rs) / 1048576))
    return out


def main():
    fixed = summarise_fixed()
    TABLES.mkdir(parents=True, exist_ok=True)
    cols = ["species", "step", "label", "n", "wall_med_h", "wall_min_h", "wall_max_h", "cores", "cpu_med_h",
            "rss_max_gb"]
    with open(TABLES / "runtime_fixed.tsv", "w") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t")
        w.writeheader()
        for r in fixed:
            w.writerow({k: f"{v:.3f}" if isinstance(v, float) else v for k, v in r.items()})
    for r in fixed:
        print(f"{r['species']:6s} {r['step']:15s} {r['label']:10s} n={r['n']:<3d} "
              f"{fmt_h(r['wall_med_h'])} ({fmt_h(r['wall_min_h'])}-{fmt_h(r['wall_max_h'])})  "
              f"CPU {r['cpu_med_h']:.2f} h  GB {r['rss_max_gb']:.1f}")
    print("wrote", TABLES / "runtime_fixed.tsv")
    rows = summarise(load())
    cols = ["species", "caller", "step", "unit", "n", "wall_med_h", "wall_min_h", "wall_max_h", "cores",
            "rss_med_gb", "rss_max_gb", "source", "note"]
    TABLES.mkdir(parents=True, exist_ok=True)
    with open(TABLES / "runtime.tsv", "w") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", restval="")
        w.writeheader()
        for r in rows:
            w.writerow({k: f"{v:.3f}" if isinstance(v, float) else v for k, v in r.items()})
    for r in rows:
        wall = (f"{fmt_h(r['wall_med_h'])} ({fmt_h(r['wall_min_h'])}-{fmt_h(r['wall_max_h'])})" if r["n"] > 1
                else fmt_h(r.get("wall_med_h", "")) or "not recorded")
        mem = "/".join(f"{r[k]:.1f}" for k in ("rss_med_gb", "rss_max_gb") if k in r) or "-"
        print(f"{r['species']:10s} {r['caller']:9s} {r['step'][:62]:62s} n={r['n']:<3d} {wall:24s} "
              f"cores {r['cores'] or '-':5s} GB {mem}")
    print("wrote", TABLES / "runtime.tsv")


if __name__ == "__main__":
    main()
