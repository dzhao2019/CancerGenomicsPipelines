# Day 5 — Quick Reference: Connecting Processes into Workflows

---

## Chaining Patterns

```groovy
// Assignment
trimmed_ch = TRIM(reads_ch)
bam_ch     = ALIGN(trimmed_ch, genome_ch)

// .out
TRIM(reads_ch)
ALIGN(TRIM.out, genome_ch)

// Pipe (single-input processes only)
Channel.fromPath("data/*.fq.gz") | COUNT_READS | view
```

---

## Named Outputs

```groovy
process FASTQC {
    output:
        tuple val(id), path("*.html"), emit: html
        tuple val(id), path("*.zip"),  emit: zip
}

FASTQC.out.html
FASTQC.out.zip
```

---

## Fan-Out / Fan-In

```groovy
// Fan-out: same channel → several processes
FASTQC(reads_ch)
TRIM(reads_ch)

// Fan-in: gather everything → one task
MULTIQC( FASTQC.out.zip.map { id, f -> f }.collect() )

// Merge streams
all_logs = FASTQC.out.zip.map { it[1] }.mix(TRIM.out.log).collect()
```

---

## Shape Matching

| Upstream emits | Downstream expects | OK? |
|---|---|---|
| `tuple val(id), path(bam)` | `tuple val(id), path(bam)` | ✅ |
| `tuple val(id), path(bam)` | `path bam` | ❌ — `.map { id, bam -> bam }` first |
| `path "*.zip"` (per task) | `path reports` (all) | ❌ — add `.collect()` |

---

## Operators Introduced

| Operator | Effect | Emits |
|---|---|---|
| `.view()` | print items | same items |
| `.collect()` | wait for all, gather into a list | 1 item (value channel) |
| `.mix(ch2)` | merge channels | all items from both |
| `.map { }` | transform each item | same count |

---

## Python ↔ Nextflow Map

| Python | Nextflow |
|---|---|
| `y = f(x)` | `y_ch = F(x_ch)` |
| `result.html` | `F.out.html` |
| `[f(x) for x in xs]` then `g(list)` | `G(F(xs_ch).collect())` |
| `a + b` | `a_ch.mix(b_ch)` |

---

## Common Mistakes

| Mistake | Symptom | Fix |
|---|---|---|
| Forgetting `.collect()` before MultiQC | MultiQC runs once per sample | `.collect()` |
| Passing a tuple to a `path` input | "Not a valid path" / staging errors | `.map` to drop the ID |
| Using `FASTQC.out` when there are several outputs | Ambiguous output error | `FASTQC.out.zip` |
| Calling the same process twice in one workflow | "Process already used" error | Alias with `include … as` (Day 15) |
| Expecting stage-by-stage execution | Confusing logs | Samples progress independently — that's the point |

---

*Day 5 · Week 1 · Nextflow Mastery Course*
