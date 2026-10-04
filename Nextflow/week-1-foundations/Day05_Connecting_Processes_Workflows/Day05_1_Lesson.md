# Day 5: Connecting Processes into Workflows
**Week 1 · Friday · 30 minutes**  
**Prerequisites**: Days 1–4 completed  
**Goal**: Chain processes into a multi-step pipeline by passing outputs as inputs

---

## 🎯 Learning Objectives

By the end of this lesson you will be able to:
1. Use a process's **output channel** as the input of the next process
2. Access outputs with assignment, `.out`, and named `emit:` outputs
3. Use the pipe operator `|` for simple linear chains
4. Fan out (one channel → several processes) and fan in (`.collect()`, `.mix()`)
5. Read a workflow block as a **dependency graph**, not a sequence of calls

---

## ⏱️ 30-Minute Breakdown

| Block | Time | Activity |
|---|---|---|
| Outputs are channels | 6 min | Chaining two processes |
| Accessing outputs | 6 min | `.out`, `emit:`, pipes |
| Fan-out / fan-in | 6 min | Reuse a channel, `collect`, `mix` |
| Hands-on exercises | 10 min | In-class exercises 1–3 |
| Reflection | 2 min | Checklist + Day 6 preview |

---

## 📖 Part 1 — Outputs Are Channels (6 min)

When a process runs, everything in its `output:` block is emitted into a **new channel**. That channel can feed the next process.

```groovy
process TRIM {
    tag "${id}"
    input:  tuple val(id), path(reads)
    output: tuple val(id), path("*_val_{1,2}.fq.gz")
    script: "trim_galore --paired ${reads}"
}

process ALIGN {
    tag "${id}"
    input:
        tuple val(id), path(reads)
        path genome
    output:
        tuple val(id), path("${id}.bam")
    script:
        """
        bwa mem ${genome} ${reads} | samtools sort -o ${id}.bam
        """
}

workflow {
    reads_ch  = Channel.fromFilePairs("data/*_R{1,2}.fastq.gz")
    genome_ch = Channel.value(file("ref/hg38.fa"))

    trimmed_ch = TRIM(reads_ch)            // output channel of TRIM
    bam_ch     = ALIGN(trimmed_ch, genome_ch)
    bam_ch.view()
}
```

```
reads_ch ─► TRIM ─► trimmed_ch ─► ALIGN ─► bam_ch
                                    ▲
                               genome_ch (value)
```

> **Key insight**: `ALIGN` for sample S1 starts **as soon as** TRIM(S1) finishes — it doesn't wait for TRIM on all samples. There are no stage barriers unless you create one.

The **output shape must match the next input shape**: TRIM emits `tuple val(id), path(...)`, ALIGN expects `tuple val(id), path(reads)` ✅.

---

## 📖 Part 2 — Accessing Outputs (6 min)

### 2.1 `.out`

```groovy
workflow {
    TRIM(reads_ch)
    ALIGN(TRIM.out, genome_ch)
}
```

### 2.2 Multiple Outputs and `emit:`

```groovy
process FASTQC {
    input:  tuple val(id), path(reads)
    output:
        tuple val(id), path("*.html"), emit: html
        tuple val(id), path("*.zip"),  emit: zip
    script: "fastqc ${reads}"
}

workflow {
    FASTQC(reads_ch)
    FASTQC.out.zip.view()      // by name ✅ (clear)
    FASTQC.out[1].view()       // by index (fragile — avoid)
}
```

### 2.3 Pipe Operator for Linear Chains

```groovy
workflow {
    Channel.fromPath("data/*.fastq.gz") | COUNT_READS | view
}
```

Pipes only work when each process has a **single** input. Once you need a reference genome, use explicit calls.

---

## 📖 Part 3 — Fan-Out and Fan-In (6 min)

### 3.1 Fan-Out: one channel, many consumers

```groovy
workflow {
    reads_ch = Channel.fromFilePairs("data/*_R{1,2}.fastq.gz")
    FASTQC(reads_ch)          // both processes receive every item
    TRIM(reads_ch)
}
```

In DSL2 a channel can be used by several processes; Nextflow copies items to each.

### 3.2 Fan-In with `.collect()` — "wait for everything"

```groovy
process MULTIQC {
    input:  path reports           // all reports together
    output: path "multiqc_report.html"
    script: "multiqc ."
}

workflow {
    FASTQC(reads_ch)
    MULTIQC( FASTQC.out.zip.map { id, zip -> zip }.collect() )
}
```

`.collect()` waits until the upstream channel is finished, then emits **one** item: a list of everything. MULTIQC runs **once**. (`.map` drops the sample ID so only files are collected — more on `map` on Day 9.)

### 3.3 Fan-In with `.mix()` — merge streams

Assume TRIM also declares `path "*_trimming_report.txt", emit: log`:

```groovy
qc_files = FASTQC.out.zip.map { it[1] }
              .mix( TRIM.out.log )
              .collect()
MULTIQC(qc_files)
```

`.mix()` interleaves items from several channels into one.

### 3.4 The Workflow Is a Graph

```
                 ┌─► FASTQC ──zip──┐
reads_ch ────────┤                 ├─mix─► collect ─► MULTIQC (1 task)
                 └─► TRIM ──log────┘
                       │
                       └─reads─► ALIGN ─► bam_ch
```

You declared the edges; Nextflow schedules the nodes.

---

## 🔗 Python ↔ Nextflow Mental Map

```python
# Python: explicit order, results held in memory
trimmed = [trim(r) for r in reads]
bams    = [align(t, genome) for t in trimmed]   # waits for ALL trimming
multiqc([fastqc(r) for r in reads])
```

```groovy
// Nextflow: declare connections; execution order is derived
workflow {
    TRIM(reads_ch)
    ALIGN(TRIM.out, genome_ch)                  // per-sample, no global barrier
    MULTIQC(FASTQC(reads_ch).zip.map{ it[1] }.collect())
}
```

| Python | Nextflow |
|---|---|
| return value | output channel |
| `y = f(x)` then `g(y)` | `G(F(x_ch))` or `F.out` |
| named tuple / dict of results | `emit:` names |
| `list(...)` then one call | `.collect()` |
| `a + b` (concatenate lists) | `.mix()` |

**Where the analogy breaks**: calling `TRIM(reads_ch)` does not run anything immediately — it adds a node to a graph. And you can call each process **only once** per workflow (Day 15 shows aliasing to reuse one).

---

## ✅ Lesson Checklist

- [ ] I can chain two processes by passing an output channel as an input
- [ ] I check that output tuple shape matches the next input
- [ ] I can use `.out`, `.out.name` and `emit:`
- [ ] I know when the pipe `|` works and when it doesn't
- [ ] I can fan out to several processes and fan in with `.collect()` / `.mix()`
- [ ] I can draw the dependency graph of a workflow block

---

## 👀 Day 6 Preview

Tomorrow: **Running Workflows and Understanding Execution**. You'll look inside `work/`, read `.command.sh` / `.command.err`, use `-resume` to skip completed tasks, and generate execution reports.
