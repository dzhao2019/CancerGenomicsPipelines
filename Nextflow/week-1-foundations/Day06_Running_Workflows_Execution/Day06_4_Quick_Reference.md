# Day 6 — Quick Reference: Running Workflows and Execution

---

## Running

```bash
nextflow run main.nf                         # run
nextflow run main.nf -resume                 # reuse cached tasks
nextflow run main.nf -resume sharp_turing    # resume a specific run
nextflow run main.nf --reads "data/*.fq.gz"  # pipeline param (double dash)
nextflow run main.nf -ansi-log false         # one line per task
nextflow run main.nf -stub-run               # use stub: blocks (no real tools)
nextflow run main.nf -w /scratch/work        # custom work dir
```

`-option` = Nextflow · `--param` = your pipeline

---

## Reports

```bash
-with-report report.html      # resource usage per process
-with-trace trace.txt          # TSV, one row per task
-with-timeline timeline.html   # Gantt chart
-with-dag flowchart.html       # pipeline graph (.html / .png / .mmd)
```

---

## Work Directory Files

| File | Contents |
|---|---|
| `.command.sh` | Your script, fully interpolated |
| `.command.run` | Wrapper: staging, env, container, then runs `.command.sh` |
| `.command.out` | stdout |
| `.command.err` | stderr |
| `.command.log` | stdout + stderr |
| `.command.trace` | CPU / memory usage |
| `.exitcode` | Exit status (0 = OK) |

```bash
cd work/a3/7f21c9*        # tab-complete the hash prefix from the console
cat .command.sh .command.err .exitcode
bash .command.run         # reproduce the task
```

---

## Logs & History

```bash
nextflow log                                       # all runs in this dir
nextflow log <run> -f name,status,exit,workdir     # fields per task
nextflow log <run> -F 'status == "FAILED"'         # filter
less .nextflow.log                                 # engine log (latest run)
```

---

## Cleaning

```bash
nextflow clean -n                 # dry run
nextflow clean -f                 # delete work dirs of the last run
nextflow clean -f -before <run>   # older runs only
nextflow clean -f -k              # remove files, keep log/metadata entries
```

---

## What Triggers Re-execution on `-resume`

| Re-runs | Cached |
|---|---|
| Process script changed | Unrelated edits elsewhere |
| Input file changed (path/size/mtime) | `publishDir` changed |
| Input value / param used in task changed | `tag` changed |
| Container changed | — |
| Any upstream task re-ran with new output | — |

`cache 'lenient'` → hash on path + size only (shared filesystems with flaky timestamps).

---

## Python ↔ Nextflow Map

| Python | Nextflow |
|---|---|
| `python run.py` | `nextflow run main.nf` |
| checkpoint files | `-resume` |
| logging module | `.nextflow.log`, `.command.err` |
| cProfile / `time` | `-with-report`, `-with-timeline` |

---

## Common Mistakes

| Mistake | Fix |
|---|---|
| Running from a different directory, then `-resume` | Launch from the same directory |
| Deleting `work/` to "save space" mid-project | Use `nextflow clean -before` |
| `-resume` typed as `--resume` | Single dash |
| No `tag` → can't tell tasks apart | `tag "${id}"` |
| Debugging from final outputs | Read `.command.sh` + `.command.err` in the task dir |

---

*Day 6 · Week 1 · Nextflow Mastery Course*
