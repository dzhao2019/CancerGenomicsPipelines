# Day 4 — In-Class Exercises: Understanding Channels
**Week 1 · Thursday · Instructor-Led**

---

## Overview

These exercises are designed for ~25 minutes of guided classroom time.  
They progress from recognition → modification → creation.  
Solutions are at the bottom — resist peeking until you've tried!

Assume this directory layout for all exercises:
```
data/
├── S1_R1.fastq.gz   S1_R2.fastq.gz
├── S2_R1.fastq.gz   S2_R2.fastq.gz
└── S3_R1.fastq.gz   S3_R2.fastq.gz
ref/
└── hg38.fa
```

---

## Exercise 1 — Predict the Output (6 min)

### Task
Write down exactly what each snippet prints (order of lines is enough).

```groovy
// A
Channel.of(1, 2, 3).view { it * 10 }

// B
Channel.fromPath("data/*_R1.fastq.gz").view { it.simpleName }

// C
Channel.fromFilePairs("data/*_R{1,2}.fastq.gz").view { id, reads -> "${id}: ${reads.size()}" }

// D
Channel.fromPath("data/*.fastq.gz").count().view()

// E
Channel.value("hg38").view()
```

<details>
<summary>💡 Solution (click to expand)</summary>

- **A**: `10`, `20`, `30`
- **B**: `S1_R1`, `S2_R1`, `S3_R1` (`simpleName` strips all extensions)
- **C**: `S1: 2`, `S2: 2`, `S3: 2`
- **D**: `6`
- **E**: `hg38`

Order between items from `fromPath`/`fromFilePairs` is not guaranteed to be alphabetical — never rely on it.
</details>

---

## Exercise 2 — Count the Tasks (6 min)

### Task
For each workflow, how many times does `ALIGN` run? Fix any that don't do what the author intended (align **every** sample).

```groovy
// A
ALIGN(Channel.fromFilePairs("data/*_R{1,2}.fastq.gz"),
      Channel.value(file("ref/hg38.fa")))

// B
ALIGN(Channel.fromFilePairs("data/*_R{1,2}.fastq.gz"),
      Channel.fromPath("ref/hg38.fa"))

// C
ALIGN(Channel.fromFilePairs("data/*_R{1,2}.fastq.gz"),
      file("ref/hg38.fa"))

// D
ALIGN(Channel.fromFilePairs("data/*_R1.fastq.gz"),
      Channel.value(file("ref/hg38.fa")))
```

<details>
<summary>💡 Solution (click to expand)</summary>

- **A**: 3 ✅ (queue of 3 + value)
- **B**: **1** ❌ — the genome is a 1-item queue channel, consumed by the first task. Fix: `Channel.value(file("ref/hg38.fa"))`.
- **C**: 3 ✅ — a plain `file(...)` passed to a process is automatically treated as a value.
- **D**: **0** ❌ — the pattern matches only R1 files, so every group holds 1 file; `fromFilePairs` expects `size: 2` by default and silently drops incomplete groups. Fix: `"data/*_R{1,2}.fastq.gz"`. Add `checkIfExists: true` so empty globs fail loudly.
</details>

---

## Exercise 3 — Build the Right Channels (7 min)

### Task
Write channel definitions for an RNA-seq pipeline needing:
1. `reads_ch` — paired-end reads as `[id, [R1, R2]]`, failing loudly if nothing matches
2. `index_ch` — a STAR index directory `ref/star_index/` reused by every sample
3. `gtf_ch` — `ref/genes.gtf` reused by every sample
4. `samples_ch` — the sample IDs `"ctrl_1", "ctrl_2", "treat_1", "treat_2"` as plain values

Add a `.view()` to `reads_ch` that prints `Sample ctrl_1 → 2 files`.

<details>
<summary>💡 Solution (click to expand)</summary>

```groovy
workflow {
    reads_ch = Channel
        .fromFilePairs("data/*_R{1,2}.fastq.gz", checkIfExists: true)
        .view { id, files -> "Sample ${id} → ${files.size()} files" }

    index_ch   = Channel.value(file("ref/star_index", type: 'dir'))
    gtf_ch     = Channel.value(file("ref/genes.gtf"))
    samples_ch = Channel.of("ctrl_1", "ctrl_2", "treat_1", "treat_2")
}
```

A directory can be passed as a single `path` input — the whole index is staged into the task.
</details>

---

## Exercise 4 — Debug: "Why did only one sample run?" (6 min)

### Task
A colleague reports that only one VCF was produced from 24 BAMs. Find the cause.

```groovy
workflow {
    bams_ch   = Channel.fromPath("bams/*.bam")
    genome_ch = Channel.fromPath("ref/hg38.fa")
    fai_ch    = Channel.fromPath("ref/hg38.fa.fai")
    dict_ch   = Channel.fromPath("ref/hg38.dict")

    HAPLOTYPECALLER(bams_ch, genome_ch, fai_ch, dict_ch)
}
```

<details>
<summary>💡 Solution (click to expand)</summary>

Three reference inputs are single-item **queue** channels. HaplotypeCaller pairs one item from each; after the first task they're empty → 1 task.

```groovy
    genome_ch = Channel.value(file("ref/hg38.fa"))
    fai_ch    = Channel.value(file("ref/hg38.fa.fai"))
    dict_ch   = Channel.value(file("ref/hg38.dict"))
```

Now HAPLOTYPECALLER runs 24 times. (On Day 10 you'll bundle FASTA + index files into one tuple.)
</details>

---

## Group Discussion Questions

1. Why do you think Nextflow made queue channels single-use instead of behaving like Python lists?
2. When might you *want* a process to run only once? Which channel type gives you that?
3. Two queue channels, one of BAMs and one of BAI index files, are passed into the same process. What could go wrong even if both have 24 items?
4. How does "a channel is a stream" relate to Nextflow's ability to start aligning sample 1 while sample 20 is still being trimmed?
