# Day 3 — Homework Exercises: Your First Nextflow Process
**Week 1 · Wednesday · Self-Paced (~45 min)**

Complete these exercises after the lesson. Solutions are at the bottom of each exercise — attempt the exercise fully before looking.

> **Test data tip**: create tiny fake FASTQs to run locally:
> ```bash
> mkdir -p data
> for s in S1 S2 S3; do
>   printf "@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nTTTTGGGG\n+\nIIIIIIII\n" | gzip > data/${s}.fastq.gz
> done
> ```

---

## Homework 1 — Wrap FastQC (15 min)

### Background
FastQC is the first step of almost every pipeline. Wrap it as a reusable process.

### Task
Create `fastqc.nf` with a process `FASTQC` that:
- Takes `tuple val(sample_id), path(reads)` as input
- Uses `tag` to show the sample ID in the log
- Uses 2 CPUs via the `cpus` directive and passes them to FastQC with `task.cpus`
- Outputs `tuple val(sample_id), path("*_fastqc.zip")` **and** `path "*_fastqc.html"`
- Publishes HTML reports to `results/fastqc` (copy mode)

Add a workflow that builds the input with:
```groovy
Channel.fromPath("data/*.fastq.gz").map { f -> [f.simpleName, f] }
```

<details>
<summary>✅ Solution</summary>

```groovy
// fastqc.nf
process FASTQC {
    tag "${sample_id}"
    cpus 2
    publishDir "results/fastqc", mode: 'copy', pattern: "*.html"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("*_fastqc.zip")
        path "*_fastqc.html"

    script:
        """
        fastqc --threads ${task.cpus} ${reads}
        """
}

workflow {
    reads_ch = Channel.fromPath("data/*.fastq.gz").map { f -> [f.simpleName, f] }
    FASTQC(reads_ch)
}
```

Notes:
- `pattern: "*.html"` restricts publishing to HTML; the ZIP stays in `work/` for downstream steps.
- `f.simpleName` strips **all** extensions (`S1.fastq.gz` → `S1`); `f.baseName` strips only the last (`S1.fastq`).
</details>

---

## Homework 2 — Fix Five Broken Processes (15 min)

### Task
Each process has one bug. Identify it and fix it.

```groovy
// A
process A_GZIP {
    input:  path txt
    output: path "${txt}.gz"
    script: "gzip -c $txt > ${txt}.gzip"
}

// B
process B_HEAD {
    input:  path fastq
    output: path "head.txt"
    script:
    '''
    zcat ${fastq} | head -n 8 > head.txt
    '''
}

// C
process C_STATS {
    input:  path bam
    output: path "${bam.simpleName}.stats"
    script:
    """
    samtools flagstat /home/alice/project/${bam} > ${bam.simpleName}.stats
    """
}

// D
process D_DATE {
    output: path "date.txt"
    script:
    """
    echo "Run at $(date)" > date.txt
    """
}

// E
process E_ALIGN {
    input:
        val sample_id
        path reads
    output:
        path "${sample_id}.bam"
    script:
    """
    bwa mem ref.fa ${reads} | samtools sort -o ${sample_id}.bam
    """
}
```

<details>
<summary>✅ Solution</summary>

- **A** — Output declares `*.gz` but script writes `*.gzip` → task fails "Missing output file". Make them match: `> ${txt}.gz`.
- **B** — `'''` single quotes: `${fastq}` is not interpolated; bash sees an empty variable. Use `"""`.
- **C** — Hard-coded absolute path bypasses staging. Use `samtools flagstat ${bam}`.
- **D** — `$(date)` inside `"""` is parsed by Groovy → compile error. Escape: `\$(date)`.
- **E** — Two problems: `ref.fa` is not declared as an input (it won't exist in the task directory), and separate `val`/`path` inputs can get **mismatched** across samples. Fix:

```groovy
process E_ALIGN {
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
(BWA also needs the index files next to `genome` — you'll handle that with multiple `path` inputs on Day 10.)
</details>

---

## Homework 3 — Python Inside a Process (15 min)

### Background
Not every step needs a command-line tool. Small Python helpers fit naturally inside processes.

### Task
Write a process `GC_CONTENT` that:
- Takes `tuple val(sample_id), path(fastq)`
- Uses a `#!/usr/bin/env python3` script block to compute the GC fraction of all reads
- Writes `${sample_id}.gc.tsv` containing `sample_id<TAB>gc_fraction`
- Emits `tuple val(sample_id), path("${sample_id}.gc.tsv")`

<details>
<summary>✅ Solution</summary>

```groovy
process GC_CONTENT {
    tag "${sample_id}"

    input:
        tuple val(sample_id), path(fastq)

    output:
        tuple val(sample_id), path("${sample_id}.gc.tsv")

    script:
    """
    #!/usr/bin/env python3
    import gzip

    gc = total = 0
    with gzip.open("${fastq}", "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seq = line.strip().upper()
                gc += seq.count("G") + seq.count("C")
                total += len(seq)

    with open("${sample_id}.gc.tsv", "w") as out:
        out.write(f"${sample_id}\\t{gc / total:.4f}\\n")
    """
}

workflow {
    Channel.fromPath("data/*.fastq.gz")
        .map { f -> [f.simpleName, f] }
        | GC_CONTENT
}
```

Watch the escaping: `${fastq}` and `${sample_id}` are Nextflow interpolations; Python's own f-string braces `{gc / total:.4f}` have no `$`, so Groovy leaves them alone; `\\t` and `\\n` produce `\t` and `\n` in the generated Python file.

Bigger Python scripts belong in the pipeline's `bin/` folder (Day 21) and are called like any tool: `gc_content.py ${fastq} > ${sample_id}.gc.tsv`.
</details>

---

## Stretch Challenge (Optional)

### Task
Add a `stub:` block to `FASTQC` from Homework 1 that just `touch`es fake outputs, then run `nextflow run fastqc.nf -stub-run`. Why is this useful when developing on a laptop without FastQC installed?

<details>
<summary>✅ Solution</summary>

```groovy
    stub:
        """
        touch ${sample_id}_fastqc.zip ${sample_id}_fastqc.html
        """
```

`-stub-run` executes the `stub:` block instead of `script:`, so you can test the **wiring** of a whole pipeline in seconds with no tools or real data. Day 25 uses this heavily for testing.
</details>

---

## Reflection Questions

Answer these in your own words (write notes in your Day 3 progress log):

1. Why does Nextflow fail a task when a declared output is missing? Why is that a feature?
2. What goes wrong if a script reads a file that is not declared as an input?
3. When would you keep a sample ID in a `tuple` rather than passing the file alone?
4. When would you put Python inside a process rather than calling a separate script?
