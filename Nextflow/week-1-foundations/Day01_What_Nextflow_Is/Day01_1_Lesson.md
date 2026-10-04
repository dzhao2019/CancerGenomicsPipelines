# Day 1: What Nextflow Actually Is
**Week 1 · Monday · 30 minutes**  
**Prerequisites**: Basic Python, familiarity with FASTQ/BAM/VCF files  
**Goal**: Understand what Nextflow is, which problem it solves, and when to choose it over a Python script

---

## 🎯 Learning Objectives

By the end of this lesson you will be able to:
1. Distinguish **data processing** (what a tool does) from **workflow orchestration** (how tools are coordinated)
2. Name the three building blocks of every Nextflow pipeline: **processes, channels, workflows**
3. Explain why Nextflow parallelises automatically while a Python `for` loop does not
4. Describe **resumability**, **portability** and **reproducibility** in one sentence each
5. Decide whether a given task is better solved with Python, Nextflow, or both

---

## ⏱️ 30-Minute Breakdown

| Block | Time | Activity |
|---|---|---|
| The coordination problem | 5 min | Why scripts break at scale |
| Core building blocks | 10 min | Processes, channels, workflows |
| The four superpowers | 5 min | Parallelism, resume, portability, reproducibility |
| Hands-on exercises | 8 min | In-class exercises 1–2 |
| Reflection | 2 min | Checklist + Day 2 preview |

---

## 📖 Part 1 — The Pipeline Coordination Problem (5 min)

A typical germline analysis for **one** sample:

```
FASTQ → FastQC → Trim Galore → BWA-MEM → samtools sort → GATK HaplotypeCaller → VCF
```

The Python version most of us write first:

```python
# pipeline.py — works fine for 3 samples
import subprocess

samples = ["S1", "S2", "S3"]
for s in samples:
    subprocess.run(f"fastqc {s}.fastq.gz", shell=True, check=True)
    subprocess.run(f"trim_galore {s}.fastq.gz", shell=True, check=True)
    subprocess.run(f"bwa mem ref.fa {s}_trimmed.fq.gz | samtools sort -o {s}.bam",
                   shell=True, check=True)
    subprocess.run(f"gatk HaplotypeCaller -R ref.fa -I {s}.bam -O {s}.vcf.gz",
                   shell=True, check=True)
```

Now scale to **500 samples on an HPC cluster**:

| Problem | What you'd have to write in Python |
|---|---|
| Samples run one after another | `multiprocessing` / job arrays / sbatch wrappers |
| Job 312 crashes at hour 30 | Checkpoint files + "skip if output exists" logic |
| Cluster uses SLURM, collaborator uses AWS | A second copy of the submission code |
| Tool versions differ between machines | Manual conda/docker management |
| "Which BAM came from which FASTQ?" | Your own naming conventions + bookkeeping |

None of this is *science* — it is **orchestration**. Nextflow is a language dedicated to orchestration so you can stop re-writing it.

> **Key insight**: Python is excellent at *processing data*. Nextflow is excellent at *coordinating the tools that process data*. Most real pipelines use both.

---

## 📖 Part 2 — The Three Building Blocks (10 min)

### 2.1 Process — one step, isolated

```groovy
process FASTQC {
    input:
        path reads

    output:
        path "*_fastqc.{html,zip}"

    script:
        """
        fastqc ${reads}
        """
}
```

A process is a **recipe**: declared inputs, declared outputs, and a shell script. Each time it runs it gets its **own private directory** under `work/`, so tasks never overwrite each other.

### 2.2 Channel — data in motion

```groovy
reads_ch = Channel.fromPath("data/*.fastq.gz")
```

A channel is a **queue of items** (files, values, tuples). It is *not* a list you loop over — items flow through it, and every item that arrives triggers a new task.

### 2.3 Workflow — the wiring diagram

```groovy
workflow {
    reads_ch = Channel.fromPath("data/*.fastq.gz")
    FASTQC(reads_ch)
}
```

The workflow block says **what connects to what**. It does **not** say "do sample 1, then sample 2". Nextflow reads the connections and schedules tasks itself.

### 2.4 Putting It Together — a 2-step pipeline

```groovy
// main.nf
process FASTQC {
    input:  path reads
    output: path "*_fastqc.zip"
    script: "fastqc ${reads}"
}

process TRIM {
    input:  path reads
    output: path "*_trimmed.fq.gz"
    script: "trim_galore ${reads}"
}

workflow {
    reads_ch = Channel.fromPath("data/*.fastq.gz")
    FASTQC(reads_ch)
    TRIM(reads_ch)
}
```

```bash
nextflow run main.nf
```

With 500 FASTQ files this launches **1000 tasks**, running as many in parallel as your machine or cluster allows — with zero loop code.

### 2.5 Declarative vs Imperative

| Imperative (Python script) | Declarative (Nextflow) |
|---|---|
| "Do A, then B, then C, for each sample" | "B needs A's output; C needs B's output" |
| You decide the execution order | Nextflow derives the order from data dependencies |
| Parallelism is extra code | Parallelism is the default |

> **The central idea of this course**: *the data-dependency graph is the program*. You describe how data flows; Nextflow decides when things run.

---

## 📖 Part 3 — The Four Superpowers (5 min)

| Superpower | What it means | Python equivalent effort |
|---|---|---|
| **Automatic parallelism** | Every channel item becomes an independent task | `multiprocessing`, job arrays |
| **Resume** (`-resume`) | Re-run skips every task whose inputs + script are unchanged | Custom checkpoint logic |
| **Portability** | Same `main.nf` on laptop, SLURM, AWS, GCP — only config changes | Rewrite submission layer |
| **Reproducibility** | Per-process containers (Docker/Singularity) pin tool versions | Manual environment management |

Quick arithmetic: 500 samples × 40 min each.
- Sequential Python: 500 × 40 = **20,000 min ≈ 14 days**
- Nextflow, 50 parallel slots: 500 / 50 × 40 = **400 min ≈ 7 h**
- Crash at sample 400, then `-resume`: only the unfinished tasks re-run.

### Where Nextflow Sits Among Alternatives

| Tool | Style | Typical home |
|---|---|---|
| **Nextflow** | Dataflow (channels push data) | nf-core, clinical genomics, HPC + cloud |
| **Snakemake** | Rule-based (pull from target files, Python-like) | Academic labs, Python-heavy groups |
| **WDL / Cromwell** | Typed task language | Broad Institute, Terra |
| **CWL** | YAML standard | Interoperability-focused projects |

---

## 🔗 Python ↔ Nextflow Mental Map

```python
# Python: imperative loop, sequential by default
for fq in glob.glob("data/*.fastq.gz"):
    run_fastqc(fq)
```

```groovy
// Nextflow: declare the connection, parallel by default
workflow {
    Channel.fromPath("data/*.fastq.gz") | FASTQC
}
```

| Python | Nextflow |
|---|---|
| function | process |
| list / generator | channel |
| `main()` calling functions in order | `workflow {}` wiring processes together |
| `subprocess.run("fastqc ...")` | `script:` block |
| checkpoint files | `-resume` |

**Where the analogy breaks**: a Python function runs the moment you call it. Calling `FASTQC(reads_ch)` in a workflow only *connects* the process to a channel — tasks start when data arrives.

**Use both together**: Nextflow orchestrates FastQC → BWA → GATK; a Python script inside a process can parse the VCF and make plots.

---

## ✅ Lesson Checklist

- [ ] I can explain orchestration vs data processing with a bioinformatics example
- [ ] I can name and describe processes, channels and workflows
- [ ] I understand that one channel item → one task → automatic parallelism
- [ ] I can explain what `-resume` saves you
- [ ] I can say when I would still use a plain Python script

---

## 👀 Day 2 Preview

Tomorrow: **Groovy Essentials for Nextflow**. Nextflow scripts are written in a Groovy-based language. You'll learn the small subset you need — string interpolation (`"${sample}"`), lists, maps and closures (`{ it * 2 }`) — mapped one-to-one onto Python syntax you already know.
