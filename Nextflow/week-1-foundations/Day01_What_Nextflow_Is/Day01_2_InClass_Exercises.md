# Day 1 — In-Class Exercises: What Nextflow Actually Is
**Week 1 · Monday · Instructor-Led**

---

## Overview

These exercises are designed for ~25 minutes of guided classroom time.  
They progress from recognition → modification → creation.  
No installation is needed today — all exercises are paper/whiteboard design work.  
Solutions are at the bottom of each exercise — resist peeking until you've tried!

---

## Exercise 1 — Python, Nextflow, or Both? (6 min)

### Task
For each scenario, decide: **Python**, **Nextflow**, or **Both** — and give a one-line reason.

| # | Scenario |
|---|---|
| A | Plot a volcano plot from one DESeq2 results table |
| B | Run FastQC + Trim Galore + STAR on 240 RNA-seq samples on the cluster |
| C | Parse one VCF and count variants per chromosome |
| D | Run BWA + GATK on 80 exomes, then summarise all VCFs into a single Excel report with pandas |
| E | Re-run a 3-day pipeline after a node failure without redoing finished samples |
| F | Prototype a new normalisation algorithm on a small count matrix |

<details>
<summary>💡 Solution (click to expand)</summary>

| # | Answer | Reason |
|---|---|---|
| A | Python | Single analysis step, no orchestration |
| B | Nextflow | Many samples × multiple tools × cluster scheduling |
| C | Python | One file, one transformation |
| D | Both | Nextflow orchestrates BWA/GATK; a Python script *inside a process* builds the report |
| E | Nextflow | `-resume` is exactly this feature |
| F | Python | Interactive prototyping; wrap it in a process later if it goes to production |

**Rule of thumb**: many samples + multiple tools + need to scale/resume → Nextflow. One file, one analysis → Python.
</details>

---

## Exercise 2 — Identify the Building Blocks (7 min)

### Task
Read the workflow below (you don't need to understand every symbol yet). Label:
1. Every **process** name
2. Every **channel**
3. Which line(s) form the **workflow wiring**
4. How many tasks run if `data/` contains 12 FASTQ files?

```groovy
process FASTQC {
    input:  path reads
    output: path "*_fastqc.zip"
    script: "fastqc ${reads}"
}

process ALIGN {
    input:
        path reads
        path genome
    output: path "*.bam"
    script: "bwa mem ${genome} ${reads} | samtools sort -o ${reads.simpleName}.bam"
}

workflow {
    reads_ch  = Channel.fromPath("data/*.fastq.gz")
    genome_ch = Channel.value(file("ref/hg38.fa"))

    FASTQC(reads_ch)
    ALIGN(reads_ch, genome_ch)
}
```

<details>
<summary>💡 Solution (click to expand)</summary>

1. Processes: `FASTQC`, `ALIGN`
2. Channels: `reads_ch` (one item per FASTQ file), `genome_ch` (a single reusable value)
3. Wiring: the four lines inside `workflow { }`
4. **24 tasks** — 12 × FASTQC + 12 × ALIGN. They all run in parallel subject to available CPUs. `genome_ch` is a *value* channel, so the same genome is reused by all 12 ALIGN tasks (you'll learn why on Day 4).
</details>

---

## Exercise 3 — Compute the Speed-up (5 min)

### Task
A WGS pipeline takes **90 min per sample**. You have **200 samples**.

1. How long does a sequential Python loop take?
2. How long does Nextflow take with **40 parallel slots** on SLURM?
3. The run fails at sample 150 (sequential) / after 150 samples are complete (Nextflow). How much work is repeated in each case if you restart?

<details>
<summary>💡 Solution (click to expand)</summary>

1. 200 × 90 = **18,000 min = 12.5 days**
2. ⌈200 / 40⌉ = 5 waves × 90 min = **450 min = 7.5 h**
3. Naïve Python restart: all 150 finished samples re-run (**225 h** wasted) unless you wrote checkpoint logic.  
   Nextflow with `-resume`: **0** finished tasks re-run; only the remaining 50 samples (2 waves ≈ 3 h).
</details>

---

## Exercise 4 — Sketch Your Own Pipeline (7 min)

### Task
Pick a pipeline you run (or want to run) in your own work — e.g. a WGS or RNA-seq pipeline on HPC.

1. List the steps (tools) in order.
2. Draw boxes for **processes** and arrows for **channels**.
3. Mark which arrows carry **one item per sample** and which carry **one shared item** (reference genome, annotation, panel of normals).
4. Mark any step that needs **all samples at once** (e.g. MultiQC, joint genotyping).

<details>
<summary>💡 Reference Solution (click to expand)</summary>

```
            ref.fa (shared)          known_sites.vcf (shared)
                 │                           │
FASTQ ──► FASTQC ┼──────────────┐            │
  │              ▼              │            ▼
  └────► TRIM ─► BWA_MEM ─► SORT ─► MARKDUP ─► BQSR ─► HAPLOTYPECALLER ─► per-sample gVCF
                                                                              │
                                              (needs ALL samples) GENOTYPE_GVCFS ◄─┘
FASTQC + MARKDUP metrics ──(needs ALL samples)──► MULTIQC
```

- Per-sample arrows: FASTQ → … → gVCF
- Shared items: `ref.fa`, `known_sites.vcf`
- All-samples steps: `GENOTYPE_GVCFS`, `MULTIQC` — on Day 5 you'll meet `.collect()`, which gathers a channel into one item for exactly this purpose.
</details>

---

## Group Discussion Questions

1. Which parts of your current scripts are **science** and which are **orchestration**?
2. What is the most painful failure you've had with a long-running script? Would `-resume` have helped?
3. Why might a clinical lab care more about **reproducibility** than raw speed?
4. Your lab already uses Snakemake. What arguments would you make for — or against — switching?
