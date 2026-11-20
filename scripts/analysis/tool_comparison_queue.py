#!/usr/bin/env python3
"""DiGAS vs SKAT comparison queue over all region files and datasets (ADNI1/2/3), through the unchanged nested pipeline.

Per region file, four arms (selection tool x protocol):
  leaky_<TOOL>_<region>   selection once on ALL subjects (paper protocol), then 10-fold classification
                          (precomputed_selection.py stands in for the selection script)
  nested_<TOOL>_<region>  selection inside each outer fold on training subjects only
TOOL = DiGAS (scripts/nested/digas3-linux_validation.py, 1000 iterations) or SKAT (scripts/analysis/skat_tool.py).

Job names carry a _ADNI2 / _ADNI3 suffix for those datasets (ADNI1 names are unsuffixed, as before).
Outputs: results/tool_comparison/<job>/ (+ <job>.log); whole-cohort selections in results/tool_comparison/selection_all/.
A job is skipped when results/tool_comparison/<job>.done exists, so the queue can be restarted after a crash.
At most --slots jobs at once (each single-threaded); a job only starts when MemAvailable > --min-free-gb.

  .venv/bin/python scripts/analysis/tool_comparison_queue.py [--slots 2] [--regions exon_5k exon ...] [--datasets ADNI1 ADNI2 ADNI3]
"""
import argparse, os, signal, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PY = str(ROOT / ".venv" / "bin" / "python")
OUT = ROOT / "results" / "tool_comparison"
SEL = OUT / "selection_all"
DATASETS = {"ADNI1": "pheno/adni1pheno.txt", "ADNI2": "pheno/adni2pheno.txt", "ADNI3": "pheno/adni3pheno.txt"}
REGIONS = ["exon_5k", "exon", "exon_protein_coding", "pre5kexon", "gene", "gene_updown20kb"]
SCRIPT = {"DiGAS": ROOT / "scripts" / "nested" / "digas3-linux_validation.py",
          "SKAT": ROOT / "scripts" / "analysis" / "skat_tool.py"}
# whole-cohort selections that already exist (exon_5k)
EXISTING = {("DiGAS", "exon_5k"): ROOT / "results" / "fold_pair_scatter_leaky" / "digas_all_subjects" / "digas_all",
            ("SKAT", "exon_5k"): ROOT / "results" / "skat_leaky" / "skat_all_exon_5k"}
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")


def mem_free_gb():
    for line in open("/proc/meminfo"):
        if line.startswith("MemAvailable"):
            return int(line.split()[1]) / 1024 ** 2
    return 0.0


def jobs_for(region, ds="ADNI1"):
    """(name, deps, argv, extra_env) in priority order."""
    J = []
    sfx = "" if ds == "ADNI1" else f"_{ds}"
    base = [PY, str(ROOT / "scripts" / "nested" / "digas_classification_nested_10fold.py"), "--root", str(ROOT),
            "--datasets", f"{ds}_QC",
            "--region-file", str(ROOT / "region" / f"{region}.txt"), "--plink", str(ROOT / "plink"), "--outer-jobs", "1"]
    for tool in ("SKAT", "DiGAS"):
        J.append((f"nested_{tool}_{region}{sfx}", [], base + ["--digas-script", str(SCRIPT[tool]),
                  "--results-dir", str(OUT / f"nested_{tool}_{region}{sfx}")], {}))
    for tool in ("SKAT", "DiGAS"):
        pre = EXISTING.get((tool, region)) if ds == "ADNI1" else None
        deps = []
        if pre is None:
            pre = SEL / f"{tool}_{region}{sfx}"
            sel = f"selection_{tool}_{region}{sfx}"
            J.append((sel, [], [PY, str(SCRIPT[tool]), "--bfile", str(ROOT / f"{ds}_QC"),
                                "--region", str(ROOT / "region" / f"{region}.txt"),
                                "--pheno", str(ROOT / DATASETS[ds]), "--out", str(pre),
                                "--iteration", "1000"], {}))
            deps = [sel]
        J.append((f"leaky_{tool}_{region}{sfx}", deps, base + ["--digas-script",
                  str(ROOT / "scripts" / "analysis" / "precomputed_selection.py"),
                  "--results-dir", str(OUT / f"leaky_{tool}_{region}{sfx}")], {"PRECOMPUTED_PREFIX": str(pre)}))
    return J


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slots", type=int, default=2)        # user limit: max 2 jobs at once (4 allowed for ADNI1-only runs)
    ap.add_argument("--max-heavy", type=int, default=2,
                    help="max concurrent DiGAS selection jobs (nested_DiGAS_*, selection_DiGAS_*: ~3 GB each)")
    ap.add_argument("--panic-gb", type=float, default=1.5,
                    help="if MemAvailable drops below this, kill the newest job started here and requeue it")
    ap.add_argument("--min-free-gb", type=float, default=4.0)
    ap.add_argument("--regions", nargs="+", default=REGIONS)
    ap.add_argument("--datasets", nargs="+", default=list(DATASETS), choices=list(DATASETS))
    ap.add_argument("--external", default="",
                    help="pgrep pattern of an outside job that occupies one slot while it runs")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); SEL.mkdir(parents=True, exist_ok=True)

    # 1) ADNI1 exon_5k, 2) ADNI1 nested SKAT of the other regions (mostly half done),
    # 3) exon_5k of the other datasets, 4) everything else, region by region, ADNI1 -> ADNI2 -> ADNI3
    alljobs = [(r, d, j) for r in a.regions for d in a.datasets for j in jobs_for(r, d)]
    tiers = [lambda r, d, n: d == "ADNI1" and r == "exon_5k",
             lambda r, d, n: d == "ADNI1" and n.startswith("nested_SKAT"),
             lambda r, d, n: r == "exon_5k",
             lambda r, d, n: True]
    order = []
    for t in tiers:
        order += [j for r, d, j in alljobs if j not in order and t(r, d, j[0])]

    done = {n for n, *_ in order if (OUT / f"{n}.done").exists()}
    skip_external = set()
    pending = [j for j in order if j[0] not in done and j[0] not in skip_external]
    if a.dry_run:
        for n, deps, argv, extra in pending:
            print(n, "deps=", deps, extra, " ".join(argv[1:3]), "...")
        return
    running, failed = {}, set()
    log = lambda m: print(f"[{time.strftime('%F %T')}] {m}", flush=True)
    heavy = lambda n: n.startswith(("nested_DiGAS", "selection_DiGAS"))

    # Adopt jobs still running from a previous queue process (killed without its children):
    # we cannot get their exit code, so success = ALL_NESTED.csv / _significant.txt written after adoption.
    adopted = {}
    for j in list(pending):
        n, deps, argv, extra = j
        key = argv[argv.index("--results-dir") + 1] if "--results-dir" in argv else argv[argv.index("--out") + 1]
        r = subprocess.run(["pgrep", "-f", "--", key + "( |$)"], stdout=subprocess.PIPE, text=True)
        pids = [int(x) for x in r.stdout.split() if int(x) != os.getpid()]
        if pids:
            pending.remove(j); adopted[n] = (min(pids), key, time.time()); log(f"ADOPT {n} (pid {min(pids)})")

    def stop_children(*_):
        for q in running.values():
            try: os.killpg(q.pid, signal.SIGTERM)
            except ProcessLookupError: pass
        sys.exit(1)
    signal.signal(signal.SIGTERM, stop_children)
    log(f"{len(pending)} jobs queued: " + ", ".join(j[0] for j in pending))

    last_mem = 0.0
    while pending or running or adopted:
        if time.time() - last_mem > 600:                        # RAM trace, to diagnose freezes
            last_mem = time.time()
            log(f"MEM free {mem_free_gb():.1f} GB, running: {', '.join(list(adopted) + list(running)) or '-'}")
        for n, (pid, key, t0) in list(adopted.items()):
            if os.path.exists(f"/proc/{pid}"):
                continue
            del adopted[n]
            outs = [Path(key) / "ALL_NESTED.csv", Path(key + "_significant.txt")]
            if any(f.exists() and f.stat().st_mtime > t0 for f in outs):
                (OUT / f"{n}.done").touch(); done.add(n); log(f"DONE {n} (adopted)")
            else:
                failed.add(n); log(f"FAILED {n} (adopted, no output), see {OUT / (n + '.log')}")
        if running and mem_free_gb() < a.panic_gb:             # emergency: free memory before the machine hangs
            n = next(reversed(running)); q = running.pop(n)
            try: os.killpg(q.pid, signal.SIGKILL)
            except ProcessLookupError: pass
            q.wait(); a.slots = max(1, a.slots - 1)
            pending.insert(0, next(j for j in order if j[0] == n))
            log(f"PANIC free RAM < {a.panic_gb} GB: killed {n}, requeued; slots now {a.slots}")
        for n, p in list(running.items()):
            rc = p.poll()
            if rc is None:
                continue
            del running[n]
            if rc == 0:
                (OUT / f"{n}.done").touch(); done.add(n); log(f"DONE {n}")
            else:
                failed.add(n); log(f"FAILED {n} (exit {rc}), see {OUT / (n + '.log')}")
        ext = bool(a.external) and subprocess.run(["pgrep", "-f", a.external], stdout=subprocess.DEVNULL).returncode == 0
        free_slots = a.slots - len(running) - len(adopted) - (1 if ext else 0)
        n_heavy = sum(heavy(n) for n in list(running) + list(adopted))
        for j in list(pending):
            if free_slots <= 0:
                break
            n, deps, argv, extra = j
            if any(d in failed for d in deps):
                pending.remove(j); failed.add(n); log(f"SKIP {n} (dependency failed)"); continue
            if not all(d in done for d in deps):
                continue
            if heavy(n) and n_heavy >= a.max_heavy:
                continue
            if mem_free_gb() < a.min_free_gb:
                break
            pending.remove(j)
            fh = open(OUT / f"{n}.log", "w")
            if n.startswith("selection_"):                     # whole-cohort selection: wall, CPU, peak RSS
                argv = ["/usr/bin/time", "-o", str(OUT / f"{n}.time"), "-f", "%e %U %S %M"] + argv
            running[n] = subprocess.Popen(argv, cwd=str(ROOT), stdout=fh, stderr=subprocess.STDOUT,
                                          env=dict(ENV1, **extra), start_new_session=True)
            free_slots -= 1; n_heavy += heavy(n)
            log(f"START {n} (free RAM {mem_free_gb():.1f} GB)")
            time.sleep(20)                                      # stagger start-up peaks
        time.sleep(10)
    log(f"queue finished; failed: {sorted(failed) or 'none'}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
