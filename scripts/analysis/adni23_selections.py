#!/usr/bin/env python3
"""Whole-cohort DiGAS and SKAT selections (exon_5k) on ADNI2 and ADNI3, run one after another as ONE extra job
next to the tool-comparison queue (user, 2026-11-16).

Same commands, names and outputs as the queue's selection_* jobs (tool_comparison_queue.jobs_for), with a .done
file each, so the queue skips them later. Guards:
  - at most 2 heavy DiGAS processes on the machine (queue's included): a DiGAS step waits while 2 are running
  - a step starts only with MemAvailable >= 6 GB (DiGAS) / 3 GB (SKAT)
  - while running, if MemAvailable < 2 GB the step is killed and retried later (after RAM recovers)

  setsid nohup .venv/bin/python scripts/analysis/adni23_selections.py >> results/tool_comparison/adni23_selections.log 2>&1 &
"""
import importlib.util, os, signal, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("Q", ROOT / "scripts/analysis/tool_comparison_queue.py")
Q = importlib.util.module_from_spec(spec); spec.loader.exec_module(Q)
ORDER = [("ADNI2", "SKAT"), ("ADNI2", "DiGAS"), ("ADNI3", "SKAT"), ("ADNI3", "DiGAS")]
log = lambda m: print(f"[{time.strftime('%F %T')}] {m}", flush=True)


def n_heavy():
    # only the python process running digas3 (not the nested script / time wrappers that mention it)
    out = subprocess.run(["pgrep", "-fc", r"^\S*python\S* \S*digas3-linux_validation\.py --bfile"],
                         capture_output=True, text=True).stdout
    return int(out.strip() or 0)


def main():
    for ds, tool in ORDER:
        name = f"selection_{tool}_exon_5k_{ds}"
        if (Q.OUT / f"{name}.done").exists():
            log(f"SKIP {name} (done)"); continue
        job = next(j for j in Q.jobs_for("exon_5k", ds) if j[0] == name)
        while True:
            need = 6.0 if tool == "DiGAS" else 3.0
            while Q.mem_free_gb() < need or (tool == "DiGAS" and n_heavy() >= 2):
                time.sleep(60)
            log(f"START {name} (free RAM {Q.mem_free_gb():.1f} GB, heavy running {n_heavy()})")
            t0 = time.time()
            with open(Q.OUT / f"{name}.log", "a") as lf:
                p = subprocess.Popen(job[2], stdout=lf, stderr=subprocess.STDOUT, cwd=ROOT, env=Q.ENV1,
                                     start_new_session=True)
                killed = False
                while p.poll() is None:
                    time.sleep(15)
                    if Q.mem_free_gb() < 2.0:
                        os.killpg(p.pid, signal.SIGKILL); p.wait(); killed = True
                        log(f"PANIC free RAM < 2 GB: killed {name}, will retry when RAM recovers")
            if killed:
                time.sleep(600); continue
            if p.returncode == 0 and Path(f"{job[2][job[2].index('--out') + 1]}_significant.txt").exists():
                (Q.OUT / f"{name}.done").touch()
                log(f"DONE {name} in {(time.time() - t0) / 60:.0f} min")
            else:
                log(f"FAILED {name} (exit {p.returncode}), see {Q.OUT / (name + '.log')}")
            break
    log("all ADNI2/ADNI3 exon_5k selections processed")


if __name__ == "__main__":
    main()
