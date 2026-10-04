# Day 5 — Homework Exercises: Connecting Processes into Workflows
**Week 1 · Friday · Self-Paced (~45 min)**

Complete these exercises after the lesson. Solutions are at the bottom of each exercise — attempt the exercise fully before looking.

> Use the fake data from Day 4 homework (`data/*_R{1,2}.fastq.gz`). Where real tools aren't installed, add `stub:` blocks (Day 3 stretch) and run with `-stub-run` to test the wiring.

---

## Homework 1 — Three-Process QC Chain (15 min)

### Task
Build `qc_chain.nf`:

1. `COUNT_READS` — input `tuple val(id), path(reads)`; writes `${id}.count.txt` containing the number of reads in R1; emits `tuple val(id), path("${id}.count.txt")`
2. `FLAG_LOW` — input from COUNT_READS; writes `${id}.flag.txt` containing `PASS` if count ≥ 1, otherwise `FAIL`
3. `SUMMARY` — input: **all** flag files; concatenates them into `summary.txt` with one line per sample: `id<TAB>PASS/FAIL`

Wire them in a workflow and print the summary file path.

<details>
<summary>✅ Solution</summary>

```groovy
// qc_chain.nf
process COUNT_READS {
    tag "${id}"
    input:  tuple val(id), path(reads)
    output: tuple val(id), path("${id}.count.txt")
    script:
        """
        echo \$(( \$(zcat ${reads[0]} | wc -l) / 4 )) > ${id}.count.txt
        """
}

process FLAG_LOW {
    tag "${id}"
    input:  tuple val(id), path(count)
    output: path "${id}.flag.txt"
    script:
        """
        n=\$(cat ${count})
        if [ "\$n" -ge 1 ]; then s=PASS; else s=FAIL; fi
        printf "${id}\\t\$s\\n" > ${id}.flag.txt
        """
}

process SUMMARY {
    input:  path flags
    output: path "summary.txt"
    script:
        """
        cat ${flags} | sort > summary.txt
        """
}

workflow {
    reads_ch = Channel.fromFilePairs("data/*_R{1,2}.fastq.gz", checkIfExists: true)
    COUNT_READS(reads_ch)
    FLAG_LOW(COUNT_READS.out)
    SUMMARY(FLAG_LOW.out.collect())
    SUMMARY.out.view { "Summary: ${it}" }
}
```

Notes: `reads[0]` picks R1 from the staged pair. Inside `"""`, bash variables are escaped (`\$n`), and `\\t` becomes a literal `\t` for `printf`.
</details>

---

## Homework 2 — Named Outputs Refactor (15 min)

### Task
Rewrite this process to use **named outputs** and update the workflow to use them. Then add a second consumer that collects only the HTML files into a process `ARCHIVE_HTML` that tars them.

```groovy
process FASTQC {
    input:  tuple val(id), path(reads)
    output:
        tuple val(id), path("*.html")
        tuple val(id), path("*.zip")
    script: "fastqc ${reads}"
}

workflow {
    FASTQC(Channel.fromFilePairs("data/*_R{1,2}.fastq.gz"))
    MULTIQC(FASTQC.out[1].map { it[1] }.collect())
}
```

<details>
<summary>✅ Solution</summary>

```groovy
process FASTQC {
    tag "${id}"
    input:  tuple val(id), path(reads)
    output:
        tuple val(id), path("*.html"), emit: html
        tuple val(id), path("*.zip"),  emit: zip
    script: "fastqc ${reads}"
}

process ARCHIVE_HTML {
    input:  path htmls
    output: path "fastqc_html.tar.gz"
    script: "tar -czf fastqc_html.tar.gz ${htmls}"
}

workflow {
    FASTQC(Channel.fromFilePairs("data/*_R{1,2}.fastq.gz"))

    MULTIQC(      FASTQC.out.zip.map  { id, f -> f }.collect() )
    ARCHIVE_HTML( FASTQC.out.html.map { id, f -> f }.collect() )
}
```

Nested lists? FastQC on a pair produces **two** HTML files per sample, so after `.map` each item is `[a.html, b.html]`. `collect()` flattens nested lists by default (`flat: true`), so the process receives one flat list of all files. Check it with `.collect().view()`.
</details>

---

## Homework 3 — Fan-Out, Fan-In Design (15 min)

### Background
A variant-calling lab wants:
- FastQC on raw reads
- Trim → BWA align → sort/index → `samtools flagstat` → GATK HaplotypeCaller
- MultiQC over FastQC zips + flagstat outputs
- `bcftools merge` over **all** VCFs into one cohort VCF

### Task
Write only the `workflow {}` block (assume all processes exist with sensible `emit:` names). Then answer: which processes run once, and which run once per sample?

<details>
<summary>✅ Solution</summary>

```groovy
workflow {
    reads_ch  = Channel.fromFilePairs(params.reads, checkIfExists: true)
    genome_ch = Channel.value(file(params.genome))

    FASTQC(reads_ch)
    TRIM(reads_ch)
    BWA_MEM(TRIM.out.reads, genome_ch)
    SAMTOOLS_INDEX(BWA_MEM.out.bam)
    FLAGSTAT(SAMTOOLS_INDEX.out.bam_bai)
    HAPLOTYPECALLER(SAMTOOLS_INDEX.out.bam_bai, genome_ch)

    MULTIQC(
        FASTQC.out.zip.map { id, f -> f }
            .mix(FLAGSTAT.out.stats.map { id, f -> f })
            .collect()
    )

    BCFTOOLS_MERGE(HAPLOTYPECALLER.out.vcf.map { id, vcf -> vcf }.collect())
}
```

- **Once per sample**: FASTQC, TRIM, BWA_MEM, SAMTOOLS_INDEX, FLAGSTAT, HAPLOTYPECALLER
- **Once**: MULTIQC, BCFTOOLS_MERGE (inputs are `.collect()`ed)

(`params.reads` / `params.genome` preview Day 8.)
</details>

---

## Stretch Challenge (Optional)

### Task
Make COUNT_READS from Homework 1 emit **two** named outputs — `count` (the file) and `value` (the read count as a value, using `env` or `stdout`). Print `S1 has 1 reads` using only the value output.

<details>
<summary>✅ Solution</summary>

```groovy
process COUNT_READS {
    tag "${id}"
    input:  tuple val(id), path(reads)
    output:
        tuple val(id), path("${id}.count.txt"), emit: count
        tuple val(id), env(N),                  emit: value
    script:
        """
        N=\$(( \$(zcat ${reads[0]} | wc -l) / 4 ))
        echo \$N > ${id}.count.txt
        """
}

workflow {
    COUNT_READS(Channel.fromFilePairs("data/*_R{1,2}.fastq.gz"))
    COUNT_READS.out.value.view { id, n -> "${id} has ${n} reads" }
}
```

`env(N)` captures a bash variable from the task as an output value — handy for small metrics.
</details>

---

## Reflection Questions

Answer these in your own words (write notes in your Day 5 progress log):

1. Why must an output's shape match the next process's input exactly?
2. When should you use `.collect()`, and what does it cost you in terms of parallelism?
3. What's the advantage of `emit:` names in a pipeline someone else will maintain?
4. How would you explain "the data-dependency graph is the program" using today's workflow?
