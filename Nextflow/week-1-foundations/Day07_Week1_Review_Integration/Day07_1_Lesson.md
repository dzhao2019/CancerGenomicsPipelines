# Day 7: Week 1 Review and Integration Project
**Week 1 · Sunday · 30 minutes**  
**Prerequisites**: Days 1–6 completed  
**Goal**: Consolidate Week 1 by reviewing every core concept and assembling them into one working multi-sample QC pipeline

---

## 🎯 Learning Objectives

By the end of this lesson you will be able to:
1. Explain the Week 1 mental model in one diagram: **processes + channels + workflow**
2. Recall the key syntax from each day without looking it up
3. Design a multi-sample QC pipeline from a written specification
4. Run, break, and resume that pipeline, and inspect its work directories
5. Identify your own weak spots before Week 2

---

## ⏱️ 30-Minute Breakdown

| Block | Time | Activity |
|---|---|---|
| Week 1 in one page | 8 min | Concept recap Days 1–6 |
| Project specification | 4 min | What you're building |
| Reference implementation walkthrough | 8 min | Reading the full pipeline |
| Hands-on | 8 min | Start in-class Exercise 1 |
| Reflection | 2 min | Self-assessment + Week 2 preview |

---

## 📖 Part 1 — Week 1 in One Page (8 min)

| Day | Concept | One-line summary | Must-know syntax |
|---|---|---|---|
| 1 | Orchestration | Nextflow coordinates tools; Python processes data | `processes · channels · workflow` |
| 2 | Groovy | Strings, lists, maps, closures | `"${x}"`, `[a:1]`, `{ it * 2 }`, `collect`, `findAll` |
| 3 | Process | Isolated unit of work with declared I/O | `input: tuple val(id), path(reads)` / `output:` / `script:` |
| 4 | Channels | Streams; queue (once) vs value (reusable) | `fromPath`, `fromFilePairs`, `Channel.value`, `.view()` |
| 5 | Workflows | Outputs are channels; wire them together | `B(A.out.x)`, `emit:`, `.collect()`, `.mix()` |
| 6 | Execution | Isolated `work/` dirs, hashing, resume, reports | `-resume`, `.command.sh`, `nextflow log`, `-with-report` |

### The Week 1 Mental Model

```
            value channel (shared reference)
                         │
FASTQ files ─► queue ─► PROCESS A ─► channel ─► PROCESS B ─► channel ─► collect ─► PROCESS C (once)
               channel   (N tasks)               (N tasks)
                 │
                 └─ each item = one task, all in parallel, each in its own work/ dir
```

> **Central insight**: you never wrote a loop this week. Parallelism came from channels, ordering came from data dependencies, and recovery came from `-resume`. *The data-dependency graph is the program.*

---

## 📖 Part 2 — Integration Project Specification (4 min)

Build `qc_pipeline.nf` that:

1. Reads **paired-end** FASTQs from `data/*_R{1,2}.fastq.gz`
2. Runs **FastQC** on every sample
3. **Counts reads** per sample and flags samples below a minimum (e.g. 1,000 reads) as `FAIL`
4. Runs **Trim Galore** only on samples that `PASS`
5. Produces **one** summary table `qc_summary.tsv` (sample, reads, status)
6. Runs **MultiQC** once over all FastQC + trimming reports
7. Is resumable and has a `tag` on every per-sample process

---

## 📖 Part 3 — Reference Implementation Walkthrough (8 min)

```groovy
// qc_pipeline.nf — Week 1 Integration Project
params.reads     = "data/*_R{1,2}.fastq.gz"
params.min_reads = 1000

// ── Processes ─────────────────────────────────────────────────────────────
process FASTQC {
    tag "${id}"
    input:  tuple val(id), path(reads)
    output: path "*_fastqc.zip", emit: zip
    script: "fastqc ${reads}"
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

process TRIM {
    tag "${id}"
    input:  tuple val(id), path(reads)
    output:
        tuple val(id), path("*_val_{1,2}.fq.gz"), emit: reads
        path "*_trimming_report.txt",             emit: log
    script: "trim_galore --paired ${reads}"
}

process SUMMARY {
    input:  val rows
    output: path "qc_summary.tsv"
    script:
        def lines = rows.collect { it.join('\t') }.join('\n')
        """
        printf "sample\\treads\\tstatus\\n" > qc_summary.tsv
        printf "${lines}\\n" >> qc_summary.tsv
        """
}

process MULTIQC {
    input:  path reports
    output: path "multiqc_report.html"
    script: "multiqc ."
}

// ── Workflow ──────────────────────────────────────────────────────────────
workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    FASTQC(reads_ch)
    COUNT_READS(reads_ch)

    // Attach a PASS/FAIL status (map = Day 9 preview)
    status_ch = COUNT_READS.out.counted.map { id, reads, n ->
        def status = n.toInteger() >= params.min_reads ? "PASS" : "FAIL"
        [id, reads, n, status]
    }

    // Trim only passing samples (filter = Day 9 preview)
    TRIM( status_ch.filter { it[3] == "PASS" }.map { id, reads, n, s -> [id, reads] } )

    SUMMARY( status_ch.map { id, reads, n, s -> [id, n, s] }.collect(flat: false) )

    MULTIQC( FASTQC.out.zip.mix(TRIM.out.log).collect() )
}
```

### What to Notice

| Line | Week 1 concept |
|---|---|
| `fromFilePairs(..., checkIfExists: true)` | Day 4 — paired-end channel factory |
| `tuple val(id), path(reads)` | Day 3 — keep ID with data |
| `env(N)` | Day 5 stretch — capture a bash value as output |
| `reads_ch` used by FASTQC and COUNT_READS | Day 5 — fan-out |
| `.collect()` before MULTIQC / SUMMARY | Day 5 — fan-in, runs once |
| `def lines = ...` above `"""` | Day 2 — Groovy code inside `script:` before the string |
| `\\t`, `\$(...)` | Day 2/3 — escaping |
| `-resume` friendly | Day 6 — every step is a pure function of its inputs |

`collect(flat: false)` keeps each `[id, n, status]` row as its own list instead of flattening everything into one long list. `params.min_reads` is a sneak peek at Day 8.

---

## 🔗 Python ↔ Nextflow Mental Map — Week 1 Summary

| Python habit | Nextflow replacement |
|---|---|
| `for sample in samples:` | A channel with one item per sample |
| `def step(x): subprocess.run(...)` | `process STEP { input / output / script }` |
| `results = step(x)` | `STEP(x_ch)` → `STEP.out` channel |
| `all_results = [..]` then summarise | `.collect()` → one task |
| `if os.path.exists(out): skip` | `-resume` |
| f-strings | GStrings `"${x}"` |
| lambda / comprehension | closure / `collect`, `findAll` |

---

## ✅ Week 1 Self-Assessment

Rate yourself 1–5 on each; anything ≤ 3 → re-read that day's Quick Reference.

- [ ] Explain why Nextflow exists (Day 1)
- [ ] Write GStrings, lists, maps and closures (Day 2)
- [ ] Write a process with tuple input and declared outputs (Day 3)
- [ ] Choose queue vs value channels correctly (Day 4)
- [ ] Chain processes, use `emit:`, `.collect()`, `.mix()` (Day 5)
- [ ] Debug via `work/` and use `-resume` (Day 6)

---

## 👀 Week 2 Preview

Week 2 turns this toy pipeline into a **practical** one: parameters and validation (Day 8), channel operators like `map`, `filter`, `groupTuple`, `join` (Day 9), multiple inputs (Day 10), `publishDir` (Day 11), error strategies (Day 12), containers (Day 13) and systematic debugging (Day 14).
