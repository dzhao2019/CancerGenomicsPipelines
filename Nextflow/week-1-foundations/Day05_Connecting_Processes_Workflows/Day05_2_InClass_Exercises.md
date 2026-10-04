# Day 5 — In-Class Exercises: Connecting Processes into Workflows
**Week 1 · Friday · Instructor-Led**

---

## Overview

These exercises are designed for ~25 minutes of guided classroom time.  
They progress from recognition → modification → creation.  
Solutions are at the bottom — resist peeking until you've tried!

---

## Exercise 1 — Draw the Graph (6 min)

### Task
Draw the dependency graph for this workflow, then answer: with 5 samples, how many tasks run in total?

```groovy
workflow {
    reads_ch  = Channel.fromFilePairs("data/*_R{1,2}.fastq.gz")
    genome_ch = Channel.value(file("ref/hg38.fa"))

    FASTQC(reads_ch)
    TRIM(reads_ch)
    ALIGN(TRIM.out.reads, genome_ch)
    FLAGSTAT(ALIGN.out)
    MULTIQC( FASTQC.out.zip.map { it[1] }
               .mix(FLAGSTAT.out.map { it[1] })
               .collect() )
}
```

<details>
<summary>💡 Solution (click to expand)</summary>

```
reads_ch ─┬─► FASTQC ─zip──────────────────┐
          │                                ├─ mix ─► collect ─► MULTIQC
          └─► TRIM ─reads─► ALIGN ─► FLAGSTAT ┘
                              ▲
                          genome_ch
```

Tasks: FASTQC 5 + TRIM 5 + ALIGN 5 + FLAGSTAT 5 + MULTIQC **1** = **21**.
</details>

---

## Exercise 2 — Fix the Shape Mismatches (7 min)

### Task
This workflow fails. Find the **three** problems and fix them.

```groovy
process TRIM {
    input:  tuple val(id), path(reads)
    output:
        tuple val(id), path("*_val_{1,2}.fq.gz"), emit: reads
        path "*_trimming_report.txt",             emit: log
    script: "trim_galore --paired ${reads}"
}

process ALIGN {
    input:
        tuple val(id), path(reads)
        path genome
    output: tuple val(id), path("${id}.bam")
    script: "bwa mem ${genome} ${reads} | samtools sort -o ${id}.bam"
}

process MULTIQC {
    input:  path reports
    output: path "multiqc_report.html"
    script: "multiqc ."
}

workflow {
    reads_ch  = Channel.fromFilePairs("data/*_R{1,2}.fastq.gz")
    genome_ch = Channel.value(file("ref/hg38.fa"))

    TRIM(reads_ch)
    ALIGN(TRIM.out, genome_ch)          // problem 1
    MULTIQC(TRIM.out.log)               // problem 2
    ALIGN.out.view { id, bam -> bam.size() } | MULTIQC   // problem 3
}
```

<details>
<summary>💡 Solution (click to expand)</summary>

1. `TRIM.out` has two outputs → ambiguous. Use `TRIM.out.reads`.
2. MultiQC would run once per sample. Use `TRIM.out.log.collect()`.
3. MULTIQC is called a second time — a process can be invoked only once per workflow. Remove that line (if you want BAM stats in MultiQC, add a FLAGSTAT process and `mix` its output into the single MULTIQC call).

```groovy
workflow {
    reads_ch  = Channel.fromFilePairs("data/*_R{1,2}.fastq.gz")
    genome_ch = Channel.value(file("ref/hg38.fa"))

    TRIM(reads_ch)
    ALIGN(TRIM.out.reads, genome_ch)
    MULTIQC(TRIM.out.log.collect())
}
```
</details>

---

## Exercise 3 — Add a Step (6 min)

### Task
Starting from the fixed workflow, add a `SAMTOOLS_INDEX` process that takes the BAM and emits `tuple val(id), path(bam), path("${bam}.bai")`. Connect it after ALIGN and `.view()` the result.

<details>
<summary>💡 Solution (click to expand)</summary>

```groovy
process SAMTOOLS_INDEX {
    tag "${id}"
    input:  tuple val(id), path(bam)
    output: tuple val(id), path(bam), path("${bam}.bai")
    script: "samtools index ${bam}"
}

workflow {
    // ... as before
    SAMTOOLS_INDEX(ALIGN.out)
    SAMTOOLS_INDEX.out.view { id, bam, bai -> "${id}: ${bam.name} + ${bai.name}" }
}
```

Re-emitting the input `bam` in the output keeps BAM and BAI travelling together — a pattern you'll use constantly.
</details>

---

## Exercise 4 — Design: Week 1 Mini Pipeline (6 min)

### Task
On paper, design the `workflow {}` block for: FastQC on raw reads → Trim Galore → FastQC on trimmed reads → MultiQC over **both** FastQC runs and the trimming logs.

Hint: you can't call `FASTQC` twice yet. Give the second one a different process name (`FASTQC_TRIMMED`) — Day 15 shows a cleaner way.

<details>
<summary>💡 Reference Solution (click to expand)</summary>

```groovy
workflow {
    reads_ch = Channel.fromFilePairs("data/*_R{1,2}.fastq.gz", checkIfExists: true)

    FASTQC(reads_ch)
    TRIM(reads_ch)
    FASTQC_TRIMMED(TRIM.out.reads)

    qc_ch = FASTQC.out.zip.map { id, f -> f }
              .mix( FASTQC_TRIMMED.out.zip.map { id, f -> f } )
              .mix( TRIM.out.log )
              .collect()

    MULTIQC(qc_ch)
}
```
</details>

---

## Group Discussion Questions

1. Why doesn't Nextflow wait for all TRIM tasks before starting any ALIGN task? When might you *want* that barrier?
2. What are the trade-offs of `.out[0]` vs `.out.reads`?
3. `.collect()` turns a queue into a value channel. Why does that matter if two processes consume it?
4. How would you explain to a Python colleague that `ALIGN(...)` inside `workflow {}` "doesn't run anything yet"?
