# Day 6 — In-Class Exercises: Running Workflows and Understanding Execution
**Week 1 · Saturday · Instructor-Led**

---

## Overview

These exercises are designed for ~25 minutes of guided classroom time.  
They progress from recognition → modification → creation.  
Solutions are at the bottom — resist peeking until you've tried!

All exercises use this pipeline (`main.nf`):

```groovy
params.reads = "data/*_R{1,2}.fastq.gz"

process COUNT {
    tag "${id}"
    input:  tuple val(id), path(reads)
    output: tuple val(id), path("${id}.n")
    script: "echo \$(( \$(zcat ${reads[0]} | wc -l) / 4 )) > ${id}.n"
}

process REPORT {
    input:  path counts
    output: path "report.tsv"
    script: "for f in ${counts}; do echo -e \"\${f%.n}\\t\$(cat \$f)\"; done > report.tsv"
}

workflow {
    COUNT(Channel.fromFilePairs(params.reads, checkIfExists: true))
    REPORT(COUNT.out.map { id, f -> f }.collect())
}
```

---

## Exercise 1 — Read the Console (5 min)

### Task
Given this output, answer the questions.

```
Launching `main.nf` [jolly_curie] DSL2 - revision: 3fa2b7c1d0
executor >  local (5)
[8d/4a1b2c] COUNT (ctrl_2)   [100%] 4 of 4 ✔
[c7/0f9e8d] REPORT           [100%] 1 of 1 ✔
```

1. What is the run name?
2. How many tasks ran in total?
3. In which directory will you find the `.command.sh` of COUNT for `ctrl_2`?
4. Can you find the work directory of COUNT for `ctrl_1` from this output alone?

<details>
<summary>💡 Solution (click to expand)</summary>

1. `jolly_curie`
2. 5 (4 COUNT + 1 REPORT)
3. `work/8d/4a1b2c…/` — the hash shown belongs to the **last-displayed** task (`ctrl_2`)
4. No — the condensed log only shows one task per process. Use `-ansi-log false` or `nextflow log jolly_curie -f tag,workdir`.
</details>

---

## Exercise 2 — Predict the Resume (8 min)

### Task
After a successful run, each change below is made **separately**, then `nextflow run main.nf -resume` is executed. For each, how many COUNT and REPORT tasks re-run?

| # | Change |
|---|---|
| A | Nothing |
| B | Edit REPORT's script (`echo -e` → `printf`) |
| C | Edit COUNT's script (divide by 4 → add a comment `# reads`) inside the script string |
| D | Add a new sample `treat_3` to `data/` |
| E | `touch data/ctrl_1_R1.fastq.gz` (update timestamp only) |
| F | Add `publishDir "results"` to REPORT |
| G | Run from a different directory with `-resume` |

<details>
<summary>💡 Solution (click to expand)</summary>

| # | COUNT re-runs | REPORT re-runs | Why |
|---|---|---|---|
| A | 0 | 0 | Everything cached |
| B | 0 | 1 | Only REPORT's script hash changed |
| C | 4 | 1 | Script text (even a comment *inside* the script) changes the hash; new outputs → REPORT re-runs |
| D | 1 | 1 | New COUNT task for `treat_3`; REPORT's collected input changed |
| E | 1 | 1 | Default caching includes timestamps → `ctrl_1` re-runs → REPORT too. (`cache 'lenient'` would avoid this) |
| F | 0 | 0 | `publishDir` isn't part of the hash — outputs are just published from cache |
| G | 4 | 1 | Different launch dir → no `.nextflow/` cache history → full re-run |

Note on D and E: with identical COUNT outputs REPORT still re-runs in D because its input list changed; in E the regenerated `ctrl_1.n` has identical content, but it's a new file in a new task directory, so REPORT's input hash changes.
</details>

---

## Exercise 3 — Investigate a Failure (7 min)

### Task
A run fails with:

```
ERROR ~ Error executing process > 'COUNT (treat_1)'
Caused by:
  Process `COUNT (treat_1)` terminated with an error exit status (1)
Command executed:
  echo $(( $(zcat treat_1_R1.fastq.gz | wc -l) / 4 )) > treat_1.n
Command exit status:
  1
Command error:
  gzip: treat_1_R1.fastq.gz: not in gzip format
Work dir:
  /home/me/proj/work/3e/91ac0d…
```

1. Which file(s) in the work dir would you inspect first, and why?
2. What is the likely root cause?
3. After fixing the input file, which command do you run, and what gets re-executed?

<details>
<summary>💡 Solution (click to expand)</summary>

1. `.command.err` (full error), `.command.sh` (exact command), and `ls -l` to see the staged `treat_1_R1.fastq.gz` symlink and where it points.
2. The input file is not actually gzipped (e.g. an uncompressed FASTQ named `.gz`, or a truncated download).
3. `nextflow run main.nf -resume` — COUNT for `treat_1` (its input changed), then REPORT. The other COUNT tasks are cached.
</details>

---

## Exercise 4 — Generate and Interpret Reports (5 min)

### Task
Write the command that resumes the pipeline and produces a report, trace and timeline. Then answer: which report would you open to find (a) the slowest task, (b) the process with the highest peak memory, (c) whether tasks actually ran in parallel?

<details>
<summary>💡 Solution (click to expand)</summary>

```bash
nextflow run main.nf -resume \
    -with-report report.html \
    -with-trace trace.txt \
    -with-timeline timeline.html
```

(a) `trace.txt` (sort by `realtime`) or the report's task table  
(b) `report.html` — memory section per process  
(c) `timeline.html` — overlapping bars = parallel tasks
</details>

---

## Group Discussion Questions

1. Why does Nextflow hash inputs + script instead of checking whether output files exist?
2. Your shared HPC filesystem changes timestamps when files are copied. What would you configure?
3. When is it safe to run `nextflow clean`? What do you lose?
4. How would `-resume` change the way you develop a pipeline compared with a Python script?
