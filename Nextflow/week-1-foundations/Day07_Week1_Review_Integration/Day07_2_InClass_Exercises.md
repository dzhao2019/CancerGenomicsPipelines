# Day 7 — In-Class Exercises: Week 1 Review and Integration Project
**Week 1 · Sunday · Instructor-Led**

---

## Overview

These exercises are designed for ~25 minutes of guided classroom time.  
Today's exercises build the Week 1 integration project step by step, then stress-test it.  
Solutions are at the bottom — resist peeking until you've tried!

> **No tools installed?** Add a `stub:` block to each process (e.g. `touch ${id}_fastqc.zip`) and run with `-stub-run`. The wiring is what matters today.

---

## Exercise 1 — Build the Skeleton (8 min)

### Task
Starting from an empty `qc_pipeline.nf`, write:
1. `params.reads` and `params.min_reads`
2. `FASTQC` (tuple input, emits `zip`)
3. `COUNT_READS` (emits `tuple val(id), path(reads), env(N)`)
4. A workflow that creates `reads_ch`, calls both processes on it, and `.view()`s COUNT_READS output

Run it and confirm you see one line per sample.

<details>
<summary>💡 Solution (click to expand)</summary>

```groovy
params.reads     = "data/*_R{1,2}.fastq.gz"
params.min_reads = 1000

process FASTQC {
    tag "${id}"
    input:  tuple val(id), path(reads)
    output: path "*_fastqc.zip", emit: zip
    script: "fastqc ${reads}"
    stub:   "touch ${id}_R1_fastqc.zip ${id}_R2_fastqc.zip"
}

process COUNT_READS {
    tag "${id}"
    input:  tuple val(id), path(reads)
    output: tuple val(id), path(reads), env(N), emit: counted
    script:
        """
        N=\$(( \$(zcat ${reads[0]} | wc -l) / 4 ))
        """
}

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    FASTQC(reads_ch)
    COUNT_READS(reads_ch)
    COUNT_READS.out.counted.view { id, reads, n -> "${id}: ${n} reads" }
}
```
</details>

---

## Exercise 2 — Add Status, Trimming and Summary (8 min)

### Task
Extend Exercise 1:
1. Add a `PASS`/`FAIL` status using `.map` (compare `n.toInteger()` with `params.min_reads`)
2. Run `TRIM` only on `PASS` samples using `.filter`
3. Add `SUMMARY` that writes `qc_summary.tsv` from `[id, n, status]` rows collected with `collect(flat: false)`
4. Add `MULTIQC` over FastQC zips + trimming logs

Test with `--min_reads 2` on the 1-read fake data — every sample should FAIL and TRIM should run **0** times.

<details>
<summary>💡 Solution (click to expand)</summary>

See the reference implementation in `Day07_1_Lesson.md`, Part 3. Key lines:

```groovy
status_ch = COUNT_READS.out.counted.map { id, reads, n ->
    [id, reads, n, n.toInteger() >= params.min_reads ? "PASS" : "FAIL"]
}
TRIM( status_ch.filter { it[3] == "PASS" }.map { id, reads, n, s -> [id, reads] } )
SUMMARY( status_ch.map { id, reads, n, s -> [id, n, s] }.collect(flat: false) )
MULTIQC( FASTQC.out.zip.mix(TRIM.out.log).collect() )
```

With every sample failing, TRIM has no input items → 0 tasks, and MULTIQC still runs on FastQC zips alone. This shows that a process simply doesn't run when its channel is empty — no `if` needed.
</details>

---

## Exercise 3 — Break, Inspect, Resume (6 min)

### Task
1. Introduce a typo into TRIM's script (`trim_galoree`). Run the pipeline.
2. From the error message, go to the failing work directory and show `.command.sh` and `.command.err`.
3. Fix the typo and resume. Which processes are cached?

<details>
<summary>💡 Solution (click to expand)</summary>

1. Error: `trim_galoree: command not found`, exit status 127.
2. `cd work/xx/yyyy…; cat .command.sh .command.err; cat .exitcode` → `127`
3. `nextflow run qc_pipeline.nf -resume` → FASTQC and COUNT_READS cached; TRIM (all samples, script changed), and MULTIQC re-run; SUMMARY cached (its inputs didn't change).

Exit code **127** always means "command not found" — usually a typo or a missing tool/container (Day 13).
</details>

---

## Exercise 4 — Code Review (3 min)

### Task
Find three improvements in this teammate's version:

```groovy
workflow {
    reads_ch = Channel.fromPath("/home/bob/data/*_R1.fastq.gz")
    FASTQC(reads_ch)
    MULTIQC(FASTQC.out.zip)
}
```

<details>
<summary>💡 Solution (click to expand)</summary>

1. Hard-coded absolute path → use `params.reads` (Day 8 makes this standard).
2. R1 only / `fromPath` → use `fromFilePairs("…_R{1,2}…")` to keep pairs and sample IDs.
3. `MULTIQC(FASTQC.out.zip)` → missing `.collect()`, so MultiQC runs once per sample.
4. (Bonus) No `checkIfExists: true` → silent empty run if the path is wrong.
</details>

---

## Group Discussion Questions

1. Which Week 1 concept was hardest, and what finally made it click?
2. In the project, where would you add a value channel if the next step were alignment?
3. Why is "a process with an empty input channel simply doesn't run" both powerful and occasionally dangerous?
4. What would you need to change to run this on 500 samples on your HPC? (Hint: most answers are Week 2 topics.)
