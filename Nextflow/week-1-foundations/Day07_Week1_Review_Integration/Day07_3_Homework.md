# Day 7 — Homework Exercises: Week 1 Review and Integration Project
**Week 1 · Sunday · Self-Paced (~45 min)**

Complete these exercises after the lesson. Solutions are at the bottom of each exercise — attempt the exercise fully before looking.

---

## Homework 1 — Finish and Harden the QC Pipeline (15 min)

### Task
Take your in-class `qc_pipeline.nf` and add:
1. A `tag` on every per-sample process
2. A `FASTQC_TRIMMED` process that runs FastQC on trimmed reads
3. MultiQC over raw FastQC + trimmed FastQC + trimming logs
4. A final `.view()` that prints the path of `qc_summary.tsv` and `multiqc_report.html`
5. Run it twice: once normally, once with `-resume`, and confirm all tasks are cached the second time

<details>
<summary>✅ Solution</summary>

```groovy
process FASTQC_TRIMMED {
    tag "${id}"
    input:  tuple val(id), path(reads)
    output: path "*_fastqc.zip", emit: zip
    script: "fastqc ${reads}"
}

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    FASTQC(reads_ch)
    COUNT_READS(reads_ch)

    status_ch = COUNT_READS.out.counted.map { id, reads, n ->
        [id, reads, n, n.toInteger() >= params.min_reads ? "PASS" : "FAIL"]
    }

    TRIM( status_ch.filter { it[3] == "PASS" }.map { id, reads, n, s -> [id, reads] } )
    FASTQC_TRIMMED(TRIM.out.reads)

    SUMMARY( status_ch.map { id, reads, n, s -> [id, n, s] }.collect(flat: false) )
    MULTIQC( FASTQC.out.zip.mix(FASTQC_TRIMMED.out.zip, TRIM.out.log).collect() )

    SUMMARY.out.view { "Summary: ${it}" }
    MULTIQC.out.view { "MultiQC: ${it}" }
}
```

Second run output should show `cached: N` for every process. `mix` accepts several channels at once.

⚠️ Raw and trimmed FastQC produce zips with **different** names (`S1_R1_fastqc.zip` vs `S1_R1_val_1_fastqc.zip`), so MultiQC staging has no name collisions. If two inputs ever share a name, Nextflow raises a "file name collision" error — Day 11 covers this.
</details>

---

## Homework 2 — Explain It Back (15 min)

### Task
Without looking at notes, write short answers (2–4 sentences each):

1. How does Nextflow decide how many FASTQC tasks to run?
2. Why is the reference genome normally a value channel?
3. What does `.collect()` do to parallelism, and why is it needed for MultiQC?
4. What exactly is in `.command.sh`, and why is it the first file you read when debugging?
5. What changes cause a task to re-run under `-resume`?

<details>
<summary>✅ Model Answers</summary>

1. One task per item in its input channel. `fromFilePairs` emits one item per sample, so FASTQC runs once per sample — no loop involved.
2. Every sample needs the same genome. A queue channel would be consumed by the first task; a value channel can be read by unlimited tasks.
3. `.collect()` waits for the upstream channel to complete and emits one list. It creates a deliberate barrier — the downstream process runs once, after everything upstream finishes. MultiQC needs every report at once.
4. The task's shell script after Nextflow interpolated all `${...}` values. It shows exactly what ran, so you can spot wrong file names, missing flags or quoting errors immediately.
5. Changes to the process script, input files (path/size/timestamp by default), input values, or container — and anything downstream of a task that re-ran with different output.
</details>

---

## Homework 3 — Plan the Week 2 Upgrade (15 min)

### Background
Your PI wants to run this QC pipeline on 300 samples on the cluster next month.

### Task
List at least **six** limitations of your current pipeline, and for each write which Week 2 day will fix it.

<details>
<summary>✅ Solution</summary>

| Limitation | Fix | Day |
|---|---|---|
| Paths and thresholds hard-coded or undocumented | `params`, validation, `--help` | 8 |
| Clumsy `it[3]` indexing; can't read a sample sheet | `map`, `filter`, `splitCsv`, `groupTuple` | 9 |
| No reference/alignment step with index files | Multiple inputs, tuples, `each` | 10 |
| Results only in `work/` | `publishDir` with organised output folders | 11 |
| One bad sample kills the whole run | `errorStrategy`, `maxRetries` | 12 |
| Requires FastQC/Trim Galore installed identically everywhere | Containers | 13 |
| Ad-hoc debugging | Systematic debugging workflow | 14 |
</details>

---

## Stretch Challenge (Optional)

### Task
Add a `-stub-run`-compatible `stub:` block to **every** process so the whole pipeline can be tested in seconds on any laptop. Then run `nextflow run qc_pipeline.nf -stub-run -with-dag dag.html` and check the graph.

<details>
<summary>✅ Solution (stubs)</summary>

```groovy
// FASTQC / FASTQC_TRIMMED
stub: "touch ${id}_R1_fastqc.zip ${id}_R2_fastqc.zip"

// COUNT_READS
stub:
    """
    N=5000
    """

// TRIM
stub:
    """
    touch ${id}_R1_val_1.fq.gz ${id}_R2_val_2.fq.gz ${id}_R1.fastq.gz_trimming_report.txt
    """

// MULTIQC
stub: "touch multiqc_report.html"
```

SUMMARY can keep its real script — it only uses `printf`. With stub N=5000 every sample PASSes, so you test the full graph.
</details>

---

## Reflection Questions

Answer these in your own words (write notes in your Week 1 section of the progress log):

1. What is the single most important idea from Week 1?
2. Which day's Quick Reference will you keep open during Week 2?
3. Which concept would you explain differently to a Python colleague now than you would have on Day 1?
4. What is one thing you want to be able to do by the end of Week 2?
