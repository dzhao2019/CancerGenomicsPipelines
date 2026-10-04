# Day 7 — Quick Reference: Week 1 Cheat Sheet

---

## Minimal Complete Pipeline

```groovy
params.reads = "data/*_R{1,2}.fastq.gz"

process FASTQC {
    tag "${id}"
    input:  tuple val(id), path(reads)
    output: path "*_fastqc.zip", emit: zip
    script: "fastqc ${reads}"
}

process MULTIQC {
    input:  path reports
    output: path "multiqc_report.html"
    script: "multiqc ."
}

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    FASTQC(reads_ch)
    MULTIQC(FASTQC.out.zip.collect())
}
```

---

## Groovy (Day 2)

```groovy
"${sample}.bam"                 // interpolation (double quotes only)
'${not_interpolated}'
[1, 2, 3].collect { it * 2 }    // [2, 4, 6]
xs.findAll { it > 1 }           // filter
[id: "S1", type: "tumor"].id    // map access
{ a, b -> a + b }               // closure with params
cond ? "yes" : "no"             // ternary
```

---

## Process (Day 3)

```groovy
process NAME {
    tag "${id}"                       // directives
    input:  tuple val(id), path(reads)
    output: tuple val(id), path("${id}.bam"), emit: bam
    script:
        """
        tool ${reads} > ${id}.bam     # \$BASH_VAR escaped
        """
    stub:   "touch ${id}.bam"
}
```

---

## Channels (Day 4)

| Factory | Emits |
|---|---|
| `Channel.of(a, b)` | a, b |
| `Channel.fromPath("*.bam")` | one Path per file |
| `Channel.fromFilePairs("*_R{1,2}.fq.gz")` | `[id, [R1, R2]]` |
| `Channel.value(file("ref.fa"))` | one reusable item |

Queue = consumed once (samples) · Value = reusable (references)

---

## Workflow (Day 5)

```groovy
B(A.out)                    // chain
B(A.out.bam)                // named output
X(ch); Y(ch)                // fan-out
Z(ch.collect())             // fan-in (runs once)
ch1.mix(ch2, ch3)           // merge
```

---

## Execution (Day 6)

```bash
nextflow run main.nf -resume
nextflow run main.nf -with-report r.html -with-trace t.txt -with-timeline tl.html
nextflow log last -f name,status,exit,workdir
cd work/ab/cdef*; cat .command.sh .command.err .exitcode
```

---

## Exit Codes Worth Knowing

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Generic tool error — read `.command.err` |
| 127 | Command not found (typo / tool or container missing) |
| 137 | Killed — usually out of memory (Day 12, 19) |
| 143 | Terminated — usually walltime exceeded / cancelled |

---

## Python ↔ Nextflow Map

| Python | Nextflow |
|---|---|
| function | process |
| list / generator | channel (queue) |
| constant | value channel |
| loop | implicit — one task per item |
| `main()` | `workflow {}` |
| checkpoints | `-resume` |
| f-string | GString |
| lambda | closure |

---

## Top 5 Week 1 Mistakes

| Mistake | Fix |
|---|---|
| Reference genome as queue channel → 1 task | `Channel.value(file(...))` |
| `'single quotes'` in `script:` | `"""triple double quotes"""` |
| Unescaped `$(cmd)` in script | `\$(cmd)` |
| Missing `.collect()` before MultiQC | add `.collect()` |
| Undeclared input file used in script | declare it with `path` |

---

*Day 7 · Week 1 · Nextflow Mastery Course*
