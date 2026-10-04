# Day 1 — Homework Exercises: What Nextflow Actually Is
**Week 1 · Monday · Self-Paced (~45 min)**

Complete these exercises after the lesson. Solutions are at the bottom of each exercise — attempt the exercise fully before looking.  
No installation is required; these are reading and design exercises.

---

## Homework 1 — Translate a Python Script into Building Blocks (15 min)

### Background
A colleague hands you this script that processes RNA-seq samples.

```python
# rnaseq.py
import glob, subprocess, os

GENOME_INDEX = "ref/star_index"
GTF          = "ref/genes.gtf"

for fq in sorted(glob.glob("data/*_R1.fastq.gz")):
    sample = os.path.basename(fq).replace("_R1.fastq.gz", "")
    r2 = fq.replace("_R1", "_R2")

    subprocess.run(f"fastqc {fq} {r2}", shell=True, check=True)
    subprocess.run(f"trim_galore --paired {fq} {r2}", shell=True, check=True)
    subprocess.run(
        f"STAR --genomeDir {GENOME_INDEX} --readFilesIn {sample}_R1_val_1.fq.gz {sample}_R2_val_2.fq.gz "
        f"--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {sample}.",
        shell=True, check=True)

subprocess.run(f"featureCounts -a {GTF} -o counts.txt *.bam", shell=True, check=True)
subprocess.run("multiqc .", shell=True, check=True)
```

### Task
1. List every **process** you would create (give each an UPPER_CASE name).
2. For each process, write its **inputs** and **outputs** in plain English.
3. Identify every **channel** and whether it carries *one item per sample* or *a single shared item*.
4. Which processes need **all samples at once**?
5. Name **three** things that would break or be painful if this script ran on 300 samples.

<details>
<summary>✅ Solution</summary>

| Process | Inputs | Outputs |
|---|---|---|
| `FASTQC` | sample's R1 + R2 | HTML/ZIP reports |
| `TRIM_GALORE` | sample's R1 + R2 | trimmed R1 + R2 |
| `STAR_ALIGN` | trimmed R1 + R2, STAR index | sorted BAM |
| `FEATURECOUNTS` | **all** BAMs, GTF | `counts.txt` |
| `MULTIQC` | **all** QC/log files | `multiqc_report.html` |

Channels:
- `reads_ch` — one item per sample (sample ID + pair of files)
- `trimmed_ch`, `bam_ch` — one item per sample
- `index_ch`, `gtf_ch` — single shared items (reused by every task)
- `all_bams_ch`, `all_qc_ch` — one item containing **all** files (gathered)

All-samples steps: `FEATURECOUNTS`, `MULTIQC`.

Painful at 300 samples:
1. Fully sequential: 300 × (QC + trim + STAR) ≈ days of wall time
2. A crash at sample 250 means re-running everything (no resume)
3. `*.bam` glob in the working directory — no isolation; stale BAMs from old runs get counted
4. (Bonus) No record of which tool versions produced the results
</details>

---

## Homework 2 — Declarative Thinking (15 min)

### Background
In Nextflow you describe **dependencies**, not **order**.

### Task
For the pipeline below, answer the questions.

```
FASTQ ─► FASTQC
FASTQ ─► TRIM ─► ALIGN ─► SORT ─► CALL_VARIANTS
                  ▲                    ▲
              genome.fa            genome.fa
```

1. With 4 samples, list one valid order in which tasks could *start* (assume unlimited CPUs).
2. Can `FASTQC` for sample 3 run at the same time as `ALIGN` for sample 1? Why?
3. Can `SORT` for sample 2 start before `TRIM` for sample 4 finishes? Why?
4. Rewrite the dependencies as a list of sentences of the form "X needs Y's output".

<details>
<summary>✅ Solution</summary>

1. Example: FASTQC×4 and TRIM×4 start immediately (8 tasks). Each ALIGN starts as soon as *its* TRIM finishes, each SORT after *its* ALIGN, etc. There is no global "stage barrier".
2. **Yes** — FASTQC(S3) depends only on FASTQ(S3); ALIGN(S1) depends only on TRIM(S1). No shared dependency.
3. **Yes** — SORT(S2) depends only on ALIGN(S2). Samples progress independently. (In a Python loop, sample 2 can't start until sample 1 is fully done.)
4. FASTQC needs FASTQ · TRIM needs FASTQ · ALIGN needs TRIM output + genome · SORT needs ALIGN output · CALL_VARIANTS needs SORT output + genome.

This list of sentences **is** your Nextflow workflow. That's the core idea: *the data-dependency graph is the program*.
</details>

---

## Homework 3 — Make the Case (15 min)

### Background
Your PI asks: "Why should we spend a month learning Nextflow instead of just improving our Python scripts?"

### Task
Write a short (≤ 200 words) answer covering:
- One concrete time saving (with numbers)
- Resumability
- Portability between your HPC and the cloud
- Reproducibility for publications / clinical reporting
- One honest limitation or cost of adopting Nextflow

<details>
<summary>✅ Example Answer</summary>

> Our current WGS scripts process samples sequentially: 100 samples × 2 h ≈ 8 days. With Nextflow on our SLURM cluster using 50 concurrent jobs the same run takes ≈ 4 h, with no parallelisation code to maintain. When a node fails — which happens on multi-day runs — `nextflow run -resume` restarts only unfinished tasks instead of everything. The same pipeline code runs on our HPC and on AWS Batch by switching a config profile, so collaborators can run it without rewriting submission scripts. Each process can pin a container image, so reviewers and auditors can reproduce exact tool versions. We also gain access to nf-core: >100 community-maintained, tested pipelines (e.g. nf-core/sarek, nf-core/rnaseq) we can use or adapt.
>
> The cost: a new language (Groovy-based DSL) and a different mental model — data flow rather than loops — which takes a few weeks to become comfortable. Debugging also moves to inspecting `work/` directories instead of print statements.
</details>

---

## Stretch Challenge (Optional)

### Task
Browse the nf-core pipeline list (https://nf-co.re/pipelines). Pick **one** pipeline related to your work (e.g. `sarek`, `rnaseq`, `methylseq`).

1. Look at its "metro map" diagram. Identify three processes and two channels.
2. Find which steps run **per sample** and which run **once for all samples**.
3. Note one feature you'd like to be able to build by Day 28.

<details>
<summary>✅ Guidance</summary>

For **nf-core/sarek**: per-sample steps include FastQC, BWA-MEM, MarkDuplicates, BQSR, HaplotypeCaller/Mutect2; once-for-all steps include MultiQC (and joint germline genotyping when enabled). Shared inputs are the reference FASTA, its indices and known-sites VCFs. You'll build a smaller RNA-seq version of exactly this structure in Week 4.
</details>

---

## Reflection Questions

Answer these in your own words (write notes in your Day 1 progress log):

1. In one sentence, what is the difference between data processing and workflow orchestration?
2. Why does declaring inputs and outputs explicitly enable both parallelism and resume?
3. Which of the four superpowers matters most for your own work, and why?
4. What part of your current workflow would you move into Nextflow first?
