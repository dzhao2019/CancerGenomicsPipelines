# Day 4 — Quick Reference: Channels

---

## Channel Factories

```groovy
Channel.of("S1", "S2", "S3")                          // explicit values
Channel.of(["S1", "tumor"], ["S2", "normal"])          // tuples
Channel.fromPath("data/*.fastq.gz")                    // files from glob
Channel.fromPath("data/*.fastq.gz", checkIfExists: true)
Channel.fromPath("data/**/*.bam")                      // recursive glob
Channel.fromFilePairs("data/*_R{1,2}.fastq.gz")        // [id, [R1, R2]]
Channel.fromFilePairs("data/*_{1,2}.fq.gz", size: 2)
Channel.value(file("ref/hg38.fa"))                     // single reusable value
Channel.empty()                                        // nothing (useful for optional inputs)
```

---

## Queue vs Value

| | Queue | Value |
|---|---|---|
| Items | many | one |
| Consumed | once each | unlimited |
| Use for | samples | references, GTF, configs |
| Created by | `of`, `fromPath`, `fromFilePairs`, process outputs | `value()`, `.collect()`, `.first()`, a plain `file()` |

---

## How Many Tasks Will Run?

| Inputs to process | Tasks |
|---|---|
| 1 queue channel with N items | N |
| queue (N) + value | N |
| queue (N) + queue (M) | min(N, M) ⚠️ |
| value + value | 1 |

---

## Inspecting Channels

```groovy
ch.view()                                   // print each item
ch.view { "Item: ${it}" }                  // custom format
ch.view { id, files -> "${id}: ${files}" } // destructure tuples
ch.count().view()                          // how many items
```

---

## Path Object Helpers

```groovy
f = file("data/S1_R1.fastq.gz")
f.name          // S1_R1.fastq.gz
f.simpleName    // S1_R1            (strips ALL extensions)
f.baseName      // S1_R1.fastq      (strips last extension)
f.extension     // gz
f.parent        // data
f.exists()      // true/false
f.size()        // bytes
```

---

## Python ↔ Nextflow Map

| Python | Nextflow |
|---|---|
| `glob.glob(pattern)` | `Channel.fromPath(pattern)` |
| generator | queue channel |
| constant | value channel |
| `print()` | `.view()` |
| `len(list)` | `.count()` |
| `for x in xs: f(x)` | `F(channel)` |

---

## Common Mistakes

| Mistake | Symptom | Fix |
|---|---|---|
| Reference as queue channel | Process runs once | `Channel.value(file(...))` |
| Glob without quotes on CLI | Shell expands it | Quote: `--reads "data/*.fq.gz"` |
| Wrong pair pattern | Samples missing / singletons | Check `{1,2}` placement; use `.view()` |
| Expecting two queue inputs to pair by sample ID | Reads matched to the wrong sample | They pair by **arrival order** — use a tuple or `join` (Day 9–10) |
| Expecting `println` order | Confusing logs | Use `.view()`; don't rely on order |

---

*Day 4 · Week 1 · Nextflow Mastery Course*
