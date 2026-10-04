# Day 3: Your First Nextflow Process
**Week 1 · Wednesday · 30 minutes**  
**Prerequisites**: Days 1–2 completed  
**Goal**: Write, run and understand a complete Nextflow process that wraps a real bioinformatics tool

---

## 🎯 Learning Objectives

By the end of this lesson you will be able to:
1. Name the parts of a process: **directives**, `input:`, `output:`, `script:`
2. Choose the right input qualifier: `val`, `path`, `tuple`
3. Declare outputs so Nextflow can find the files a tool produced
4. Run a single-process workflow and find its results in `work/`
5. Explain why a process is *not* a Python function

---

## ⏱️ 30-Minute Breakdown

| Block | Time | Activity |
|---|---|---|
| Process anatomy | 5 min | The four parts of a process |
| Inputs and outputs | 8 min | `val`, `path`, `tuple`; output globs |
| Script block & running | 5 min | Interpolation, escaping, `nextflow run` |
| Hands-on exercises | 10 min | In-class exercises 1–3 |
| Reflection | 2 min | Checklist + Day 4 preview |

---

## 📖 Part 1 — Process Anatomy (5 min)

```groovy
process FASTQC {                              // ① name (UPPER_CASE by convention)
    tag "${reads.simpleName}"                 // ② directives (optional settings)
    publishDir "results/fastqc", mode: 'copy'

    input:                                    // ③ what the task receives
        path reads

    output:                                   // ④ what the task produces
        path "*_fastqc.{html,zip}"

    script:                                   // ⑤ the shell command(s)
        """
        fastqc ${reads}
        """
}
```

| Part | Purpose | Python analogy |
|---|---|---|
| Name | How the workflow calls it | function name |
| Directives | Settings: `tag`, `publishDir`, `cpus`, `memory`, `container` | decorators |
| `input:` | Declared inputs, staged into the task directory | parameters |
| `output:` | Files/values collected after the script finishes | return value |
| `script:` | Bash run in an isolated directory | function body (`subprocess.run`) |

---

## 📖 Part 2 — Inputs and Outputs (8 min)

### 2.1 Input Qualifiers

```groovy
input:
    val  sample_id                    // a plain value (string, number, map)
    path reads                        // a file — Nextflow stages (symlinks) it into the task dir
    tuple val(sample_id), path(reads) // several related items kept together
```

> **Why `path` matters**: each task runs in its own directory `work/ab/12cd…/`. Declaring `path reads` tells Nextflow to put a link to that file there. A file you reference by absolute path in the script without declaring it is invisible to resume, containers and cloud executors.

### 2.2 Output Declarations

```groovy
output:
    path "${sample_id}.bam"                       // exact name
    path "*.{html,zip}"                           // glob
    tuple val(sample_id), path("*.vcf.gz")        // keep the ID attached to the file
    path "stats.txt", emit: stats                 // named output (used on Day 5)
```

Nextflow checks for the declared files after `script:` finishes. **If a declared output is missing, the task fails** — a built-in sanity check Python scripts don't give you.

### 2.3 Tuple Pattern — keep sample ID and data together

```groovy
process ALIGN {
    tag "${sample_id}"

    input:
        tuple val(sample_id), path(reads)
        path genome

    output:
        tuple val(sample_id), path("${sample_id}.bam")

    script:
        """
        bwa mem ${genome} ${reads} | samtools sort -o ${sample_id}.bam
        """
}
```

This `tuple val(id), path(file)` pattern is how every serious pipeline tracks "which BAM came from which sample". In Week 3 `val(sample_id)` grows into `val(meta)` — a map with id, condition, etc.

---

## 📖 Part 3 — The Script Block and Running a Process (5 min)

### 3.1 Interpolation and Escaping (Day 2 recap)

```groovy
script:
    """
    echo "Sample ${sample_id} on \$(hostname)"   # Nextflow var / escaped bash
    fastqc --threads ${task.cpus} ${reads}
    """
```

`task.cpus` is the value of the `cpus` directive — use it instead of hard-coding thread counts.

### 3.2 A Complete Runnable Script

```groovy
// main.nf
process COUNT_READS {
    tag "${fastq.simpleName}"

    input:
        path fastq

    output:
        path "${fastq.simpleName}.count.txt"

    script:
        """
        echo \$(( \$(zcat ${fastq} | wc -l) / 4 )) > ${fastq.simpleName}.count.txt
        """
}

workflow {
    COUNT_READS(Channel.fromPath("data/*.fastq.gz"))
}
```

```bash
nextflow run main.nf
```

Expected output:
```
executor >  local (3)
[3f/a1b2c3] COUNT_READS (S1) [100%] 3 of 3 ✔
```

`[3f/a1b2c3]` is the start of the task's work directory: `work/3f/a1b2c3…/`. Inside you'll find `S1.count.txt` plus hidden files `.command.sh`, `.command.out`, `.command.err`, `.exitcode` (Day 6 explores these).

### 3.3 Other Script Forms

```groovy
script:
    """
    #!/usr/bin/env python3
    import gzip
    n = sum(1 for _ in gzip.open("${fastq}")) // 4
    open("${fastq.simpleName}.count.txt", "w").write(str(n))
    """
```

A shebang line lets you write the task body in **Python** (or R). Nextflow still handles staging, parallelism and resume.

---

## 🔗 Python ↔ Nextflow Mental Map

```python
# Python
def count_reads(fastq: str) -> str:
    out = fastq.replace(".fastq.gz", ".count.txt")
    n = int(subprocess.check_output(f"zcat {fastq} | wc -l", shell=True)) // 4
    open(out, "w").write(str(n))
    return out
```

```groovy
// Nextflow
process COUNT_READS {
    input:  path fastq
    output: path "${fastq.simpleName}.count.txt"
    script:
    """
    echo \$(( \$(zcat ${fastq} | wc -l) / 4 )) > ${fastq.simpleName}.count.txt
    """
}
```

| Python function | Nextflow process |
|---|---|
| Runs when called | Runs once **per channel item** when data arrives |
| Shares the current directory | Gets its own isolated `work/` directory |
| Returns any object | Outputs must be declared |
| Silent if a file isn't written | Task **fails** if a declared output is missing |
| Called many times anywhere | Can be *invoked* once per workflow scope (aliasing on Day 15) |

---

## ✅ Lesson Checklist

- [ ] I can label name, directives, input, output and script in any process
- [ ] I know when to use `val`, `path` and `tuple`
- [ ] I understand why files must be declared with `path`, not hard-coded
- [ ] I can escape bash variables inside the script block
- [ ] I ran (or traced) a one-process workflow and know where `work/` outputs live
- [ ] I can explain three differences between a process and a Python function

---

## 👀 Day 4 Preview

Tomorrow: **Understanding Channels**. Today you fed a process with `Channel.fromPath(...)` without looking closely. Tomorrow you'll learn queue vs value channels, `Channel.of`, `fromFilePairs`, and why channels — not loops — are what make Nextflow parallel.
