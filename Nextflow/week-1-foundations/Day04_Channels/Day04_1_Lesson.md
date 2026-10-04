# Day 4: Understanding Channels
**Week 1 · Thursday · 30 minutes**  
**Prerequisites**: Days 1–3 completed  
**Goal**: Understand channels as streams of data, create them with channel factories, and see how they drive automatic parallelism

---

## 🎯 Learning Objectives

By the end of this lesson you will be able to:
1. Explain a channel as a **stream** (conveyor belt), not a list
2. Distinguish **queue channels** (consumed once) from **value channels** (reused forever)
3. Create channels with `Channel.of`, `Channel.fromPath`, `Channel.fromFilePairs`, `Channel.value`
4. Inspect channel contents with `.view()`
5. Predict how many tasks a process will run given its input channels

---

## ⏱️ 30-Minute Breakdown

| Block | Time | Activity |
|---|---|---|
| What is a channel? | 4 min | Stream vs list |
| Queue vs value | 6 min | Consumption rules, the classic "only one task ran" bug |
| Channel factories | 8 min | `of`, `fromPath`, `fromFilePairs`, `value` |
| Hands-on exercises | 10 min | In-class exercises 1–3 |
| Reflection | 2 min | Checklist + Day 5 preview |

---

## 📖 Part 1 — What Is a Channel? (4 min)

```python
# Python: load everything, then loop sequentially
files = glob.glob("data/*.fastq.gz")
for f in files:
    run_fastqc(f)
```

```groovy
// Nextflow: items flow in; each item triggers a task
workflow {
    Channel.fromPath("data/*.fastq.gz") | FASTQC
}
```

Think of a channel as a **conveyor belt**. Items are placed on it; a process downstream picks each item up as it arrives and starts a task. Nothing waits for "the whole list".

```
fromPath ──[S1.fq]──[S2.fq]──[S3.fq]──►  FASTQC  ──► 3 parallel tasks
```

> **Key insight**: The number of tasks a process runs is determined by **how many items arrive on its input channel(s)** — not by any loop you write.

---

## 📖 Part 2 — Queue vs Value Channels (6 min)

| | Queue channel | Value channel |
|---|---|---|
| Created by | `Channel.of`, `fromPath`, `fromFilePairs`, process outputs | `Channel.value(x)`, `file(x)` passed directly, `.collect()`, `.first()` |
| Holds | Many items, in order | Exactly one item |
| Read | Each item consumed **once** | Can be read **unlimited** times |
| Typical use | Samples | Reference genome, GTF, config files |

### The Classic Bug: "Only one sample was aligned"

```groovy
workflow {
    reads_ch  = Channel.fromPath("data/*.fastq.gz")       // 10 items (queue)
    genome_ch = Channel.fromPath("ref/hg38.fa")           // 1 item  (queue!) ❌
    ALIGN(reads_ch, genome_ch)
}
// → ALIGN runs ONCE. After the first task consumed the genome, the queue is empty.
```

When a process has several queue inputs, it pairs items **one-from-each** and stops when the shortest channel runs out.

```groovy
workflow {
    reads_ch  = Channel.fromPath("data/*.fastq.gz")       // 10 items
    genome_ch = Channel.value(file("ref/hg38.fa"))        // value ✅
    ALIGN(reads_ch, genome_ch)
}
// → ALIGN runs 10 times; every task gets the same genome.
```

> Rule: **per-sample data → queue channel; shared reference data → value channel.**

---

## 📖 Part 3 — Channel Factories (8 min)

### 3.1 `Channel.of` — explicit values

```groovy
Channel.of("S1", "S2", "S3").view()
// S1
// S2
// S3

Channel.of(["S1", "tumor"], ["S2", "normal"]).view()
// [S1, tumor]
// [S2, normal]
```

### 3.2 `Channel.fromPath` — files from a glob

```groovy
Channel.fromPath("data/*.fastq.gz").view()
// /abs/path/data/S1.fastq.gz
// /abs/path/data/S2.fastq.gz

Channel.fromPath("data/*.fastq.gz", checkIfExists: true)   // fail if no match
```

Items are **Path objects**, so you can call `.simpleName`, `.baseName`, `.extension`, `.size()`.

### 3.3 `Channel.fromFilePairs` — paired-end reads

```groovy
Channel.fromFilePairs("data/*_R{1,2}.fastq.gz").view()
// [S1, [/abs/data/S1_R1.fastq.gz, /abs/data/S1_R2.fastq.gz]]
// [S2, [/abs/data/S2_R1.fastq.gz, /abs/data/S2_R2.fastq.gz]]
```

Each item is already a `[sample_id, [R1, R2]]` tuple — exactly what `tuple val(sample_id), path(reads)` expects (Day 3).

### 3.4 `Channel.value` — a single reusable item

```groovy
genome_ch = Channel.value(file("ref/hg38.fa"))
gtf_ch    = Channel.value(file("ref/genes.gtf"))
```

### 3.5 `.view()` — your print statement

```groovy
Channel.fromFilePairs("data/*_R{1,2}.fastq.gz")
    .view { id, files -> "Sample ${id} has ${files.size()} files" }
```

`.view()` passes items through unchanged, so you can insert it anywhere in a chain while debugging.

### 3.6 Channels Are Asynchronous

```groovy
Channel.of(1, 2, 3).view()
println "Hello"
// Output order may be: Hello, 1, 2, 3
```

Channel operations are **declarations** executed by the dataflow engine; top-level `println` runs immediately. Don't rely on print ordering.

> 📝 Newer Nextflow versions also accept lowercase `channel.of(...)`. Both spellings work; this course uses `Channel.`.

---

## 🔗 Python ↔ Nextflow Mental Map

| Python | Nextflow |
|---|---|
| `glob.glob("*.fq.gz")` | `Channel.fromPath("*.fq.gz")` |
| list of tuples `[(id, [r1, r2]), ...]` | `Channel.fromFilePairs("*_R{1,2}.fq.gz")` |
| generator (consumed once) | queue channel |
| constant / global variable | value channel |
| `print(x)` | `.view()` |
| `for f in files: process(f)` | `PROCESS(channel)` — one task per item, in parallel |

**Where the analogy breaks**: a Python list can be iterated many times; a queue channel cannot. And a channel isn't "filled first, then used" — producers and consumers run at the same time.

---

## ✅ Lesson Checklist

- [ ] I can explain why a channel is a stream, not a list
- [ ] I know the difference between queue and value channels
- [ ] I can explain (and fix) the "only one task ran" bug
- [ ] I can create channels with `of`, `fromPath`, `fromFilePairs`, `value`
- [ ] I use `.view()` to inspect channel contents
- [ ] I can predict the number of tasks from channel sizes

---

## 👀 Day 5 Preview

Tomorrow: **Connecting Processes into Workflows**. A process's **output** is itself a channel. You'll chain FASTQC → TRIM → ALIGN by passing outputs as inputs, use `.out` and named `emit:` outputs, and use `.collect()` to gather all samples for MultiQC.
