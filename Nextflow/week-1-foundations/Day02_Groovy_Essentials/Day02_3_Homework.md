# Day 2 — Homework Exercises: Groovy Essentials
**Week 1 · Tuesday · Self-Paced (~45 min)**

Complete these exercises after the lesson. Solutions are at the bottom of each exercise — attempt the exercise fully before looking.

> **How to run Groovy snippets**: save them as `test.nf` and run `nextflow run test.nf`. A Nextflow script with no processes simply executes the Groovy code top to bottom. (`groovy test.groovy` also works if Groovy is installed.)

---

## Homework 1 — Sample Sheet Wrangling (15 min)

### Background
You received a sample list as Groovy data (later in the course this will come from a CSV).

```groovy
def rows = [
    [sample: "P01_T", patient: "P01", type: "tumor",  depth: 82],
    [sample: "P01_N", patient: "P01", type: "normal", depth: 35],
    [sample: "P02_T", patient: "P02", type: "tumor",  depth: 28],
    [sample: "P02_N", patient: "P02", type: "normal", depth: 31],
    [sample: "P03_T", patient: "P03", type: "tumor",  depth: 95],
]
```

### Task
Write Groovy code that prints:
1. The list of all sample names
2. Only tumour samples with depth ≥ 30
3. A list of BAM file names in the form `P01_T.sorted.bam`
4. The mean depth across all samples (rounded to 1 decimal)
5. A map from patient → list of sample names, e.g. `[P01:[P01_T, P01_N], ...]`

<details>
<summary>✅ Solution</summary>

```groovy
// 1
def names = rows.collect { it.sample }
println names                                  // [P01_T, P01_N, P02_T, P02_N, P03_T]

// 2
def good_tumours = rows.findAll { it.type == "tumor" && it.depth >= 30 }
println good_tumours.collect { it.sample }     // [P01_T, P03_T]

// 3
def bams = rows.collect { "${it.sample}.sorted.bam" }
println bams

// 4
def mean = rows.collect { it.depth }.sum() / rows.size()
println String.format("%.1f", mean as double)  // 54.2

// 5
def by_patient = rows.groupBy { it.patient }
                     .collectEntries { patient, rs -> [patient, rs.collect { it.sample }] }
println by_patient   // [P01:[P01_T, P01_N], P02:[P02_T, P02_N], P03:[P03_T]]
```

**Python equivalents**: (1) list comprehension, (2) comprehension with `if`, (4) `statistics.mean`, (5) `itertools.groupby` / `defaultdict(list)`. Groovy's `groupBy` + `collectEntries` replaces the `defaultdict` loop.
</details>

---

## Homework 2 — Fix the Quoting Bugs (15 min)

### Background
Most Groovy bugs in Nextflow are quoting bugs. The script block below should print the sample name, the current date (from bash), and write a file named after the sample.

### Task
Find and fix **four** bugs.

```groovy
def sample_id = "NA12878"
def threads   = 8

def script = '''
echo "Sample: ${sample_id}"
echo "Date: $(date)"
bwa mem -t $threads ref.fa ${sample_id}.fq.gz > ${sample_id}.sam
echo "Done with ${sample_ID}"
'''
println script
```

<details>
<summary>✅ Solution</summary>

1. `'''` (single-quoted triple string) does **not** interpolate → use `"""`.
2. Once you switch to `"""`, `$(date)` is parsed as a Groovy expression → escape it: `\$(date)`.
3. `$threads` works in Groovy, but `${threads}` is clearer and safer when followed by other characters — use braces consistently.
4. `${sample_ID}` — Groovy is case-sensitive; variable is `sample_id`.

```groovy
def sample_id = "NA12878"
def threads   = 8

def script = """
echo "Sample: ${sample_id}"
echo "Date: \$(date)"
bwa mem -t ${threads} ref.fa ${sample_id}.fq.gz > ${sample_id}.sam
echo "Done with ${sample_id}"
"""
println script
```

**Rule**: inside `"""…"""`, `${...}` belongs to Nextflow/Groovy; anything meant for bash gets a backslash: `\$VAR`, `\$(cmd)`.
</details>

---

## Homework 3 — Parse File Names into Tuples (15 min)

### Background
Tomorrow's processes will take inputs like `[sample_id, file]`. Building these tuples from file names is pure Groovy.

### Task
Given:
```groovy
def files = [
    "data/NA12878_S1_L001_R1_001.fastq.gz",
    "data/NA12878_S1_L001_R2_001.fastq.gz",
    "data/HG002_S2_L001_R1_001.fastq.gz",
    "data/HG002_S2_L001_R2_001.fastq.gz",
]
```

Write code that produces:
```
[[NA12878, [data/NA12878_S1_L001_R1_001.fastq.gz, data/NA12878_S1_L001_R2_001.fastq.gz]],
 [HG002,   [data/HG002_S2_L001_R1_001.fastq.gz,   data/HG002_S2_L001_R2_001.fastq.gz]]]
```

Hints: extract the base name after the last `/`; the sample ID is the first `_`-separated token; group by ID; sort each group so R1 comes before R2.

<details>
<summary>✅ Solution</summary>

```groovy
def pairs = files
    .groupBy { f -> f.tokenize("/")[-1].tokenize("_")[0] }   // id → [files]
    .collect { id, fs -> [id, fs.sort()] }                     // [[id, [R1, R2]], ...]

pairs.each { println it }
```

In Nextflow you'll rarely do this by hand: `Channel.fromFilePairs("data/*_R{1,2}_001.fastq.gz")` (Day 4) does exactly this for you. Knowing the Groovy underneath makes the channel factory far less mysterious.
</details>

---

## Stretch Challenge (Optional)

### Task
Write a closure `make_meta` that takes a file name like `P01_tumor_rep2.fastq.gz` and returns a meta map:
```groovy
[id: "P01_tumor_rep2", patient: "P01", condition: "tumor", replicate: 2]
```
Then apply it to a list of five file names with `collect` and print only the `tumor` entries.

<details>
<summary>✅ Solution</summary>

```groovy
def make_meta = { String fname ->
    def base  = fname.replaceAll(/\.fastq\.gz$/, "")
    def parts = base.tokenize("_")
    [id: base, patient: parts[0], condition: parts[1],
     replicate: parts[2].replace("rep", "") as Integer]
}

def names = ["P01_tumor_rep1.fastq.gz", "P01_tumor_rep2.fastq.gz",
             "P01_normal_rep1.fastq.gz", "P02_tumor_rep1.fastq.gz",
             "P02_normal_rep1.fastq.gz"]

names.collect(make_meta).findAll { it.condition == "tumor" }.each { println it }
```

Note `/regex/` slashy strings — Groovy's raw-string syntax for regular expressions (like Python `r"..."`).
</details>

---

## Reflection Questions

Answer these in your own words (write notes in your Day 2 progress log):

1. Why does the difference between `'...'` and `"..."` matter so much more in Nextflow than in Python?
2. When would you use `it` vs a named closure parameter? Which is more readable in a pipeline other people maintain?
3. Which Groovy collection method feels most different from Python, and how will you remember it?
4. Look at `.map { file -> [file.simpleName, file] }`. Can you now explain every symbol in it?
