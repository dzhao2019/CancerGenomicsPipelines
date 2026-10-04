# Day 2: Groovy Essentials for Nextflow
**Week 1 · Tuesday · 30 minutes**  
**Prerequisites**: Day 1 completed  
**Goal**: Learn just enough Groovy to read and write Nextflow scripts — mapped directly onto Python you already know

---

## 🎯 Learning Objectives

By the end of this lesson you will be able to:
1. Use **string interpolation** (`"${sample}.bam"`) and know when single quotes disable it
2. Create and manipulate **lists** and **maps** (Python lists and dicts)
3. Write **closures** (`{ x -> x * 2 }`, `{ it * 2 }`) — Groovy's lambdas
4. Use `collect`, `findAll`, `each` the way you use comprehensions and `filter` in Python
5. Recognise Groovy's optional parentheses, which make Nextflow code look "magic"

---

## ⏱️ 30-Minute Breakdown

| Block | Time | Activity |
|---|---|---|
| Why Groovy? | 2 min | Where Groovy appears in a Nextflow script |
| Strings & variables | 6 min | Interpolation, quoting, multi-line strings |
| Lists, maps, closures | 10 min | Core data structures and lambdas |
| Hands-on exercises | 10 min | In-class exercises 1–3 |
| Reflection | 2 min | Checklist + Day 3 preview |

---

## 📖 Part 1 — Why Groovy? (2 min)

Nextflow's language is built on **Groovy** (a JVM language). You only need a small subset:

```groovy
params.outdir = "results"                        // variable assignment

Channel.fromPath("data/*.fastq.gz")              // method call
    .map { file -> [file.simpleName, file] }     // closure (lambda)
    .view()

script:
"""
bwa mem ${genome} ${reads} > ${sample_id}.sam    // string interpolation
"""
```

> **Key insight**: 90% of the Groovy you'll write is string interpolation, lists, maps and closures. Skip classes, inheritance, and the rest of the Groovy book.

---

## 📖 Part 2 — Variables and Strings (6 min)

### 2.1 Variables

```groovy
def sample  = "NA12878"      // def = local variable (like Python assignment)
def depth   = 30             // Integer
def ratio   = 0.95           // BigDecimal (decimal)
def paired  = true           // Boolean
def nothing = null           // Python None
```

Semicolons are optional; types are inferred.

### 2.2 String Interpolation — the most important pattern

```groovy
def sample = "NA12878"
def reads  = "${sample}_R1.fastq.gz"     // "NA12878_R1.fastq.gz"
def msg    = "Depth is ${30 * 2}x"       // expressions allowed
def short_ = "Sample: $sample"           // braces optional for simple names
```

| Python | Groovy |
|---|---|
| `f"{sample}.bam"` | `"${sample}.bam"` |
| `'literal {x}'` | `'literal ${x}'` — single quotes do **not** interpolate |
| `"""multi\nline"""` | `"""multi\nline"""` — triple double quotes **do** interpolate |

> ⚠️ **Nextflow consequence**: process `script:` blocks use `"""…"""`, so `${reads}` is filled in by Nextflow. A **bash** variable must be escaped: `\$HOME`, `\$(date)`.

```groovy
script:
"""
echo "Processing ${sample_id}"          # Nextflow variable → substituted
echo "Started at \$(date) on \$HOSTNAME"  # bash → escaped with backslash
"""
```

### 2.3 Useful String Methods

```groovy
def f = "NA12878_R1.fastq.gz"
f.replace("_R1", "_R2")       // "NA12878_R2.fastq.gz"
f.tokenize("_")[0]            // "NA12878"      (like split, drops empties)
f.endsWith(".gz")             // true
f.toUpperCase()               // "NA12878_R1.FASTQ.GZ"
```

---

## 📖 Part 3 — Lists, Maps and Closures (10 min)

### 3.1 Lists (Python lists)

```groovy
def samples = ["S1", "S2", "S3"]
samples[0]              // "S1"
samples[-1]             // "S3"
samples.size()          // 3   (len(samples))
samples << "S4"         // append
samples.contains("S2")  // true   ("S2" in samples)
```

### 3.2 Maps (Python dicts)

```groovy
def meta = [id: "S1", condition: "tumor", single_end: false]
meta.id                 // "S1"   (meta["id"] also works)
meta.condition = "normal"
meta.keySet()           // [id, condition, single_end]
def bigger = meta + [batch: 2]   // merge → new map
```

> Keys are strings by default — `[id: "S1"]` is the same as Python `{"id": "S1"}`. This **meta map** pattern becomes central in Weeks 3–4.

### 3.3 Closures (lambdas)

```groovy
def double_it = { x -> x * 2 }      // explicit parameter
double_it(5)                         // 10

def triple_it = { it * 3 }          // implicit parameter "it"
triple_it(5)                         // 15

def label = { id, cond -> "${id}_${cond}" }   // two parameters
label("S1", "tumor")                          // "S1_tumor"
```

### 3.4 Collection Methods (comprehensions)

| Python | Groovy |
|---|---|
| `[x*2 for x in xs]` | `xs.collect { it * 2 }` |
| `[x for x in xs if x > 10]` | `xs.findAll { it > 10 }` |
| `next(x for x in xs if x > 10)` | `xs.find { it > 10 }` |
| `for x in xs: print(x)` | `xs.each { println it }` |
| `any(x > 10 for x in xs)` | `xs.any { it > 10 }` |
| `sum(xs)` | `xs.sum()` |
| `", ".join(xs)` | `xs.join(", ")` |

```groovy
def files = ["S1_R1.fq.gz", "S2_R1.fq.gz", "S3_R1.fq.gz"]
def ids   = files.collect { it.tokenize("_")[0] }        // [S1, S2, S3]
def pairs = files.collect { f -> [f.tokenize("_")[0], f] }
// [[S1, S1_R1.fq.gz], [S2, S2_R1.fq.gz], [S3, S3_R1.fq.gz]]
```

### 3.5 Optional Parentheses

```groovy
println("hello")      // same as
println "hello"

reads_ch.map({ it.simpleName })   // same as
reads_ch.map { it.simpleName }    // closure as last argument goes outside ()
```

That is why Nextflow code like `publishDir "results", mode: 'copy'` is just a method call with arguments.

### 3.6 Control Flow (unchanged from Python ideas)

```groovy
if (depth >= 30) { println "OK" } else { println "Low coverage" }

def status = depth >= 30 ? "pass" : "fail"     // ternary
def name   = meta.alias ?: meta.id             // Elvis: "alias or id"
def len    = meta?.id?.size()                  // safe navigation (no NPE on null)
```

---

## 🔗 Python ↔ Nextflow Mental Map

```python
# Python
samples = ["S1", "S2", "S3"]
bams = [f"{s}.sorted.bam" for s in samples if s != "S2"]
meta = {"id": "S1", "condition": "tumor"}
print(f"{meta['id']} is {meta['condition']}")
```

```groovy
// Groovy
def samples = ["S1", "S2", "S3"]
def bams = samples.findAll { it != "S2" }.collect { "${it}.sorted.bam" }
def meta = [id: "S1", condition: "tumor"]
println "${meta.id} is ${meta.condition}"
```

**Same concepts**: variables, f-strings, lists, dicts, lambdas, comprehensions.  
**Where the analogy breaks**: single vs double quotes matter; `it` is an implicit parameter with no Python equivalent; maps use `[:]`, not `{}` (braces mean a closure).

---

## ✅ Lesson Checklist

- [ ] I can interpolate variables with `"${x}"` and know single quotes don't interpolate
- [ ] I know to escape bash variables inside `"""` script blocks (`\$VAR`)
- [ ] I can create and index lists and maps
- [ ] I can write closures with an explicit parameter and with `it`
- [ ] I can translate a Python list comprehension into `collect` / `findAll`
- [ ] I recognise optional parentheses and trailing closures

---

## 👀 Day 3 Preview

Tomorrow: **Your First Nextflow Process**. You'll combine today's string interpolation with `input:`, `output:` and `script:` blocks to wrap a real tool (FastQC) as a portable, isolated unit of work.
