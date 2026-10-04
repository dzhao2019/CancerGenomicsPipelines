# Day 6 — Homework Exercises: Running Workflows and Understanding Execution
**Week 1 · Saturday · Self-Paced (~45 min)**

Complete these exercises after the lesson. Solutions are at the bottom of each exercise — attempt the exercise fully before looking.

> Use the `main.nf` from today's in-class exercises and the fake data from Day 4.

---

## Homework 1 — Break It, Fix It, Resume It (15 min)

### Task
1. Add a process `SLOW` between COUNT and REPORT that sleeps 10 s, then copies its input to `${id}.slow.n`. Make it **fail** for sample `treat_2` only (hint: `if [ "${id}" == "treat_2" ]; then exit 1; fi`).
2. Run the pipeline and confirm it fails. Record: how many SLOW tasks completed?
3. Remove the failure line and run with `-resume`. Record which tasks were cached.
4. Use `nextflow log` to print `name,status,exit,duration` for both runs.

<details>
<summary>✅ Solution</summary>

```groovy
process SLOW {
    tag "${id}"
    input:  tuple val(id), path(n)
    output: tuple val(id), path("${id}.slow.n")
    script:
        """
        sleep 10
        if [ "${id}" == "treat_2" ]; then exit 1; fi
        cp ${n} ${id}.slow.n
        """
}

workflow {
    COUNT(Channel.fromFilePairs(params.reads, checkIfExists: true))
    SLOW(COUNT.out)
    REPORT(SLOW.out.map { id, f -> f }.collect())
}
```

1st run: COUNT 4/4 ✔; SLOW fails on `treat_2`. By default (`errorStrategy 'terminate'`) Nextflow kills running tasks when the first failure occurs, so 0–3 other SLOW tasks may have completed depending on timing. REPORT never runs.

2nd run (`-resume`, failure line removed): **all** SLOW tasks re-run, because you edited SLOW's script — the hash changed for every sample. COUNT is cached (4/4).

> Lesson: to keep successful tasks cached, fix *data* or *upstream* problems without editing the process script — or accept the re-run. Day 12 shows `errorStrategy 'finish'` / `'ignore'` to let other samples complete.

```bash
nextflow log                                   # get the two run names
nextflow log <run1> -f name,status,exit,duration
nextflow log <run2> -f name,status,exit,duration
```
</details>

---

## Homework 2 — Work Directory Scavenger Hunt (15 min)

### Task
For the REPORT task of your last successful run:
1. Find its work directory **without** scrolling the console (use `nextflow log`).
2. Show that its inputs are symlinks and print where one points.
3. Show the interpolated script. What did `${counts}` become?
4. Re-run the task by hand with `bash .command.run` and confirm `report.tsv` is regenerated.
5. Explain in one sentence why re-running by hand is such a powerful debugging tool.

<details>
<summary>✅ Solution</summary>

```bash
# 1
nextflow log last -f name,workdir | grep REPORT
cd <that workdir>

# 2
ls -l                     # *.n files shown as -> /…/work/xx/…/ctrl_1.n
readlink -f ctrl_1.slow.n

# 3
cat .command.sh
# for f in ctrl_1.slow.n ctrl_2.slow.n treat_1.slow.n treat_2.slow.n; do ...
# ${counts} became a space-separated list of the staged file names

# 4
rm report.tsv
bash .command.run
cat report.tsv

# 5
```
Re-running `.command.run` reproduces the task in exactly the same staged environment (inputs, container, env), so you can iterate on a failing command without re-running the pipeline.
</details>

---

## Homework 3 — Execution Reports (15 min)

### Task
1. Run the pipeline from scratch (no `-resume`) with report, trace, timeline and DAG.
2. From `trace.txt`, which task had the longest `realtime`?
3. From `timeline.html`, did the four SLOW tasks overlap? What limits how many run at once on your machine?
4. Open `flowchart.html`. Does the graph match the one you would draw by hand?
5. Configure tracing in a `nextflow.config` so you never have to type the flags (preview of Day 20).

<details>
<summary>✅ Solution</summary>

```bash
nextflow run main.nf \
    -with-report report.html -with-trace trace.txt \
    -with-timeline timeline.html -with-dag flowchart.html

sort -t$'\t' -k<realtime column> trace.txt   # or open in a spreadsheet
```

2. One of the SLOW tasks (~10 s).
3. Yes — they overlap if you have ≥ 4 CPUs. The local executor runs as many tasks as fit in available CPUs (each task requests 1 CPU by default); `executor.cpus` or `maxForks` can limit it.
4. COUNT → SLOW → collect → REPORT.
5.
```groovy
// nextflow.config
trace    { enabled = true; file = "pipeline_info/trace.txt"; overwrite = true }
report   { enabled = true; file = "pipeline_info/report.html"; overwrite = true }
timeline { enabled = true; file = "pipeline_info/timeline.html"; overwrite = true }
dag      { enabled = true; file = "pipeline_info/dag.html"; overwrite = true }
```
</details>

---

## Stretch Challenge (Optional)

### Task
Your HPC copies input data nightly, updating timestamps and invalidating `-resume`. Add one directive to COUNT to make its caching robust to this, and explain the risk of doing so.

<details>
<summary>✅ Solution</summary>

```groovy
process COUNT {
    cache 'lenient'
    ...
}
```

`lenient` hashes input files by **path and size** only (ignoring timestamps). Risk: if a file's content changes but its size stays identical, the stale cached result is reused. `cache 'deep'` hashes file *content* — safest but slow for large files.
</details>

---

## Reflection Questions

Answer these in your own words (write notes in your Day 6 progress log):

1. What are the three files you will check first in any failed task directory?
2. Why does editing a process's script re-run all of its tasks on `-resume`?
3. Which execution report would you show a PI who asks "why does this pipeline take 3 days"?
4. What habits will you adopt to keep `-resume` working (launch directory, `work/`, inputs)?
