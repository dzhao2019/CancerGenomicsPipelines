# Day 6: Running Workflows and Understanding Execution
**Week 1 · Saturday · 30 minutes**  
**Prerequisites**: Days 1–5 completed  
**Goal**: Run pipelines confidently, understand the `work/` directory, and use `-resume` and execution reports

---

## 🎯 Learning Objectives

By the end of this lesson you will be able to:
1. Run a pipeline with `nextflow run` and read the live console output
2. Navigate a task's work directory and explain each hidden `.command.*` file
3. Use `-resume` and predict which tasks will be cached
4. Use `nextflow log` to find past runs and failed tasks
5. Generate `-with-report`, `-with-trace`, `-with-timeline` and `-with-dag` outputs

---

## ⏱️ 30-Minute Breakdown

| Block | Time | Activity |
|---|---|---|
| Running & console output | 5 min | `nextflow run`, reading the progress lines |
| The work directory | 7 min | Task hashes, `.command.*` files |
| Resume | 6 min | How caching works, what invalidates it |
| Hands-on exercises | 10 min | In-class exercises 1–3 |
| Reflection | 2 min | Checklist + Day 7 preview |

---

## 📖 Part 1 — Running a Pipeline (5 min)

```bash
nextflow run main.nf                       # local script
nextflow run main.nf --reads "data/*.fq.gz"   # pass a param (Day 8)
nextflow run nf-core/rnaseq -profile test,docker --outdir out   # remote pipeline (Day 16)
```

> Single dash `-resume` = **Nextflow option**. Double dash `--reads` = **pipeline parameter**.

Console output:

```
N E X T F L O W  ~  version 24.10.0
Launching `main.nf` [sharp_turing] DSL2 - revision: 1a2b3c4d5e

executor >  local (9)
[a3/7f21c9] FASTQC (S1)          [100%] 3 of 3 ✔
[5b/0e44d1] TRIM (S3)            [100%] 3 of 3 ✔
[e1/9c3a20] ALIGN (S2)           [ 66%] 2 of 3
```

| Piece | Meaning |
|---|---|
| `sharp_turing` | Random **run name** (used by `nextflow log`, `-resume <name>`) |
| `executor > local (9)` | 9 tasks submitted to the local executor so far |
| `[a3/7f21c9]` | Work-directory prefix of the **last** task shown for that process |
| `(S1)` | The `tag` directive value — always set one! |
| `3 of 3 ✔` | Completed / total tasks for that process |

Add `-ansi-log false` to print one line per task instead of the condensed view.

---

## 📖 Part 2 — The Work Directory (7 min)

Every task runs in its own directory named by a **hash** of its inputs + script:

```
work/
└── a3/
    └── 7f21c9e0b8d24f...      ← one task (FASTQC on S1)
        ├── S1_R1.fastq.gz -> /abs/data/S1_R1.fastq.gz   (staged input: symlink)
        ├── S1_R1_fastqc.zip                              (output)
        ├── .command.sh       ← the exact script after interpolation
        ├── .command.run      ← wrapper Nextflow submits (staging, env, container)
        ├── .command.out      ← stdout
        ├── .command.err      ← stderr
        ├── .command.log      ← stdout + stderr combined
        ├── .command.trace    ← resource usage (cpu, memory)
        └── .exitcode         ← 0 = success
```

### Debugging Recipe (expanded on Day 14)

```bash
cd work/a3/7f21c9e0b8d24f*
cat .command.sh          # what actually ran?
cat .command.err         # what went wrong?
cat .exitcode            # how did it exit?
bash .command.run        # re-run the task exactly as Nextflow did
```

> **Key insight**: `.command.sh` shows your script with every `${...}` already substituted. Most "why didn't this work" questions are answered by reading it.

---

## 📖 Part 3 — Resume (6 min)

```bash
nextflow run main.nf            # first run — ALIGN fails on S3
# fix the bug...
nextflow run main.nf -resume    # FASTQC, TRIM, ALIGN(S1,S2) cached; only ALIGN(S3)+ re-run
```

```
[a3/7f21c9] FASTQC (S1)  [100%] 3 of 3, cached: 3 ✔
[5b/0e44d1] TRIM (S3)    [100%] 3 of 3, cached: 3 ✔
[f0/12ab34] ALIGN (S3)   [100%] 3 of 3, cached: 2 ✔
```

### What Invalidates the Cache?

A task's hash includes:

| Changes that **re-run** a task | Changes that do **not** |
|---|---|
| Script text of the process | Comments/whitespace elsewhere in the file |
| Input file content/path/size/timestamp | Other processes' scripts (unless their outputs change) |
| Input values (`val`, params used in the script) | Output directory (`publishDir`) |
| Container image / conda env | Directives like `tag` |

Changing an upstream task re-runs it **and everything downstream of it**.

### Resume Gotchas

- Resume needs the **same launch directory** (it reads `.nextflow/` cache there) and the `work/` directory.
- Deleting `work/` = no resume.
- Input files touched/modified (e.g. re-copied) get new timestamps → cache miss. `cache 'lenient'` hashes path + size only.
- `-resume <run_name>` resumes from a specific earlier run.

---

## 📖 Part 4 — Logs and Reports (bonus reading)

```bash
nextflow log                                  # list past runs
nextflow log sharp_turing -f name,status,exit,workdir
nextflow log sharp_turing -F 'status == "FAILED"'   # only failed tasks
cat .nextflow.log                             # detailed engine log of the last run

# report = resources per process · trace = one TSV line per task
# timeline = Gantt chart of tasks · dag = pipeline graph
nextflow run main.nf -resume \
    -with-report report.html \
    -with-trace trace.txt \
    -with-timeline timeline.html \
    -with-dag flowchart.html

nextflow clean -n                             # preview deleting old work dirs
nextflow clean -f -before sharp_turing        # delete work dirs of older runs
```

---

## 🔗 Python ↔ Nextflow Mental Map

| Python script | Nextflow |
|---|---|
| `python pipeline.py` | `nextflow run main.nf` |
| Outputs in current dir (may overwrite) | Each task isolated in `work/xx/hash` |
| `print` / logging | `.command.out`, `.command.err`, `.nextflow.log` |
| Re-run from scratch after a crash | `-resume` |
| `time python ...` / profiler | `-with-report`, `-with-timeline`, `-with-trace` |
| "Skip if output exists" hacks | Content-based task hashing |

**Where the analogy breaks**: Python checkpoint logic usually checks "does the output file exist?". Nextflow checks "are the *inputs and script* identical?" — so a stale output from changed code is never reused by mistake.

---

## ✅ Lesson Checklist

- [ ] I can read the console progress lines (run name, hash, tag, counts)
- [ ] I can find a task's work directory and read `.command.sh` / `.command.err`
- [ ] I can predict what `-resume` will re-run after a change
- [ ] I know what breaks resume (moving launch dir, deleting `work/`, touching inputs)
- [ ] I can list past runs and failed tasks with `nextflow log`
- [ ] I can generate report, trace, timeline and DAG files

---

## 👀 Day 7 Preview

Tomorrow: **Week 1 Review and Integration Project**. You'll combine everything — Groovy, processes, channels, workflows and execution — into a complete multi-sample QC pipeline, run it, break it, and resume it.
