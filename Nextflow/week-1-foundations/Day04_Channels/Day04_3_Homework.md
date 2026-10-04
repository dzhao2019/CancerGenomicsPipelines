# Day 4 — Homework Exercises: Understanding Channels
**Week 1 · Thursday · Self-Paced (~45 min)**

Complete these exercises after the lesson. Solutions are at the bottom of each exercise — attempt the exercise fully before looking.

> **Setup**: create fake paired-end data:
> ```bash
> mkdir -p data ref
> for s in ctrl_1 ctrl_2 treat_1 treat_2; do
>   for r in 1 2; do
>     printf "@r\nACGT\n+\nIIII\n" | gzip > data/${s}_R${r}.fastq.gz
>   done
> done
> echo ">chr1" > ref/genome.fa; echo "ACGTACGTACGT" >> ref/genome.fa
> ```

---

## Homework 1 — Channel Explorer Script (15 min)

### Task
Create `explore.nf` (no processes) that prints:
1. Every R1 file's `simpleName`
2. The total number of FASTQ files
3. Each pair as `ctrl_1 → ctrl_1_R1.fastq.gz, ctrl_1_R2.fastq.gz` (file **names** only, not full paths)
4. The genome file name from a value channel, printed **twice** by viewing the same channel twice

Run with `nextflow run explore.nf`.

<details>
<summary>✅ Solution</summary>

```groovy
// explore.nf
workflow {
    // 1
    Channel.fromPath("data/*_R1.fastq.gz")
        .view { "R1: ${it.simpleName}" }

    // 2
    Channel.fromPath("data/*.fastq.gz")
        .count()
        .view { "Total FASTQ files: ${it}" }

    // 3
    Channel.fromFilePairs("data/*_R{1,2}.fastq.gz", checkIfExists: true)
        .view { id, reads -> "${id} → ${reads*.name.join(', ')}" }

    // 4
    genome_ch = Channel.value(file("ref/genome.fa"))
    genome_ch.view { "Genome (1st read): ${it.name}" }
    genome_ch.view { "Genome (2nd read): ${it.name}" }
}
```

`reads*.name` is Groovy's **spread operator**: it calls `.name` on every element (like `[r.name for r in reads]`).
</details>

---

## Homework 2 — Queue vs Value Experiment (15 min)

### Background
Seeing the "only one task ran" bug yourself makes it stick.

### Task
Create `queue_vs_value.nf` with a process `PAIR_WITH_REF` that takes `tuple val(id), path(reads)` and `path ref`, and outputs `stdout` with `echo "${id} aligned to ${ref}"`.

Run it **twice**:
1. With `ref` as `Channel.fromPath("ref/genome.fa")`
2. With `ref` as `Channel.value(file("ref/genome.fa"))`

Record how many tasks ran each time and explain the difference in 2–3 sentences.

<details>
<summary>✅ Solution</summary>

```groovy
// queue_vs_value.nf
params.mode = "value"   // try: --mode queue

process PAIR_WITH_REF {
    tag "${id}"
    input:
        tuple val(id), path(reads)
        path ref
    output:
        stdout
    script:
        """
        echo "${id} aligned to ${ref}"
        """
}

workflow {
    reads_ch = Channel.fromFilePairs("data/*_R{1,2}.fastq.gz")
    ref_ch   = params.mode == "queue"
               ? Channel.fromPath("ref/genome.fa")
               : Channel.value(file("ref/genome.fa"))

    PAIR_WITH_REF(reads_ch, ref_ch).view()
}
```

```bash
nextflow run queue_vs_value.nf --mode queue   # 1 task
nextflow run queue_vs_value.nf --mode value   # 4 tasks
```

**Explanation**: A process with multiple queue inputs takes one item from each per task. `fromPath` produced a queue with a single item, which the first task consumed, so no further tasks could be formed. A value channel is never consumed, so it pairs with every read pair.

(`params.mode` is a sneak peek at Day 8; `output: stdout` captures what the script prints.)
</details>

---

## Homework 3 — Choose the Factory (15 min)

### Task
For each situation, write the channel definition you'd use and say whether it's queue or value.

1. Single-end FASTQs: `raw/*.fq.gz`
2. Paired FASTQs named `SampleA_1.fq.gz` / `SampleA_2.fq.gz`
3. A BED file of target regions used by every sample
4. A list of chromosomes to call variants on separately: chr1–chr22, chrX
5. All BAM files in nested folders `runs/*/bams/*.bam`
6. A known-sites VCF and its `.tbi` index, used by every sample

<details>
<summary>✅ Solution</summary>

```groovy
// 1 — queue
Channel.fromPath("raw/*.fq.gz", checkIfExists: true)

// 2 — queue
Channel.fromFilePairs("raw/*_{1,2}.fq.gz", checkIfExists: true)

// 3 — value
Channel.value(file("ref/targets.bed"))

// 4 — queue (one task per chromosome later)
Channel.of(*(1..22).collect { "chr${it}" }, "chrX")

// 5 — queue
Channel.fromPath("runs/*/bams/*.bam")

// 6 — value (a single item that is a list of two files)
Channel.value([file("ref/dbsnp.vcf.gz"), file("ref/dbsnp.vcf.gz.tbi")])
```

Note on 4: `*` spreads a list into separate arguments, so `Channel.of` emits 23 items instead of 1 list. On Day 10 you'll combine this chromosome channel with BAMs to parallelise variant calling by region.
</details>

---

## Stretch Challenge (Optional)

### Task
Write a channel that emits `[sample_id, condition, [R1, R2]]` for the four samples, where `condition` is `ctrl` or `treat` derived from the sample ID. Use `fromFilePairs` plus a `.map {}` closure (Day 2 Groovy + a preview of Day 9).

<details>
<summary>✅ Solution</summary>

```groovy
Channel.fromFilePairs("data/*_R{1,2}.fastq.gz")
    .map { id, reads -> [id, id.tokenize("_")[0], reads] }
    .view()
// [ctrl_1, ctrl, [.../ctrl_1_R1.fastq.gz, .../ctrl_1_R2.fastq.gz]]
// ...
```
</details>

---

## Reflection Questions

Answer these in your own words (write notes in your Day 4 progress log):

1. Describe a queue channel and a value channel using a non-computing analogy.
2. Why is `checkIfExists: true` worth adding to almost every input channel?
3. What would you check first if a process ran fewer tasks than expected?
4. How does treating data as a stream enable parallelism without a loop?
