# Day 4: Channels and Workflows - Hands-On Exercises

**Time**: 20 minutes for exercises + 15 minutes for solutions  
**Difficulty**: Easy to Moderate  
**Goal**: Build and understand Nextflow workflows

---

## Exercise 1: Understanding Channels (5 minutes)

### Part A: Create Channels

**What does each of these create?**

```groovy
// 1.
ch1 = Channel.fromPath("data/*.fastq")

// 2.
ch2 = Channel.from("sample1", "sample2", "sample3")

// 3.
ch3 = Channel.value("reference.fa")

// 4.
ch4 = Channel.fromFilePairs("reads/*_{1,2}.fastq.gz")
```

**For each, describe:**
- What type of channel is it?
- What does it emit?
- How many times would a process run if given this channel?

<details>
<summary>✓ Solution A</summary>

**1. `Channel.fromPath("data/*.fastq")`**
- **Type:** Queue channel
- **What it emits:** File paths matching the pattern
- **Process runs:** Once for each matching file (if 5 .fastq files, process runs 5 times)
- **Example:** If files are sample1.fastq, sample2.fastq, sample3.fastq
  - Emits 3 items: /path/sample1.fastq, /path/sample2.fastq, /path/sample3.fastq

**2. `Channel.from("sample1", "sample2", "sample3")`**
- **Type:** Queue channel
- **What it emits:** The three string values
- **Process runs:** 3 times (once per string)
- **Example:** Emits "sample1", "sample2", "sample3" sequentially

**3. `Channel.value("reference.fa")`**
- **Type:** Value channel
- **What it emits:** Single value (reference file path)
- **Process runs:** Once per item in other channel (used by all)
- **Example:** Reference is available to all processes that need it

**4. `Channel.fromFilePairs("reads/*_{1,2}.fastq.gz")`**
- **Type:** Queue channel with tuples
- **What it emits:** Tuples of [sample_id, [file1, file2]]
- **Process runs:** Once per pair (if 3 paired samples, runs 3 times)
- **Example:** 
  - [sample1, [reads/sample1_1.fastq.gz, reads/sample1_2.fastq.gz]]
  - [sample2, [reads/sample2_1.fastq.gz, reads/sample2_2.fastq.gz]]

**Key insight:** Queue channels parallelize (each item separately). Value channels broadcast (same value to all).

</details>

---

### Part B: Choose the Right Channel

**For each scenario, which channel should you use?**

| Scenario | Channel Type |
|----------|--------------|
| You have 100 FASTQ files to process | Channel.from... or Channel.fromPath? |
| You need a reference genome for all samples | Channel.value... or Channel.fromPath? |
| You have paired-end reads (R1 and R2) | Channel.fromFilePairs or Channel.fromPath? |
| You have a list of sample names | Channel.from or Channel.fromPath? |

<details>
<summary>✓ Solution B</summary>

| Scenario | Choice | Why |
|----------|--------|-----|
| 100 FASTQ files | `Channel.fromPath` | Discovers files dynamically, parallelizes each |
| Reference for all | `Channel.value` | Same file used by all processes |
| Paired-end reads | `Channel.fromFilePairs` | Keeps pairs together as units |
| Sample name list | `Channel.from` | Direct list values |

</details>

---

## Exercise 2: Transform Channels (5 minutes)

### Part A: Using .map()

**What does this channel transformation produce?**

```groovy
files = Channel.fromPath("data/*.fastq")
names = files.map { file -> file.baseName }
```

**If `data/` contains:**
- `tumor_001.fastq`
- `normal_001.fastq`
- `tumor_002.fastq`

**What does `names` channel emit?**

<details>
<summary>✓ Solution A</summary>

**The names channel emits:**
```
tumor_001
normal_001
tumor_002
```

**Explanation:**
- `.map { file -> ... }` transforms each item
- `file.baseName` returns filename without extension
- So "data/tumor_001.fastq" becomes "tumor_001"

**More complex example:**
```groovy
tuples = files.map { file ->
    [file.baseName, file]  // Returns [name, path] tuple
}

// tuples emits:
// [tumor_001, /data/tumor_001.fastq]
// [normal_001, /data/normal_001.fastq]
// [tumor_002, /data/tumor_002.fastq]
```

</details>

---

### Part B: Using .filter()

**What does this produce?**

```groovy
files = Channel.fromPath("data/*.fastq")
large_files = files.filter { file -> file.size() > 1000000 }
```

**If files have these sizes:**
- tumor_001.fastq: 500 KB
- normal_001.fastq: 2 MB
- tumor_002.fastq: 1.5 MB

**What does `large_files` emit?**

<details>
<summary>✓ Solution B</summary>

**The large_files channel emits:**
```
/data/normal_001.fastq    (2 MB > 1 MB)
/data/tumor_002.fastq     (1.5 MB > 1 MB)
```

**Explanation:**
- `.filter { ... }` keeps only items where condition is true
- `file.size() > 1000000` checks if file is larger than 1 MB
- tumor_001.fastq (500 KB) is filtered out

**Note:** Sizes are in bytes!
- 1 MB = 1,000,000 bytes
- 1 GB = 1,000,000,000 bytes

</details>

---

## Exercise 3: Building a Simple Workflow (10 minutes)

### Your Task

**Write a workflow that:**
1. Finds all FASTQ files in `data/`
2. Creates sample names from filenames
3. Passes to a FASTQC process (already defined)
4. Views the results

**Your workflow:**

```groovy
workflow {
    ???
}
```

<details>
<summary>✓ Solution</summary>

```groovy
workflow {
    // Step 1: Find all FASTQ files
    samples = Channel.fromPath("data/*.fastq")
    
    // Step 2: Transform to sample names and tuples
    named_samples = samples.map { file ->
        [file.baseName, file]
    }
    
    // Step 3: Pass to process
    qc_results = FASTQC(named_samples)
    
    // Step 4: View results
    qc_results.view()
}
```

**What happens:**
1. `Channel.fromPath` finds all FASTQ files
2. `.map` creates tuples like [sample_id, file_path]
3. `FASTQC` process receives each tuple
4. Runs once per sample (in parallel)
5. Results viewed at the end

**Alternative (more concise):**
```groovy
workflow {
    samples = Channel.fromPath("data/*.fastq")
        .map { f -> [f.baseName, f] }
    
    FASTQC(samples).view()
}
```

**Even more realistic:**
```groovy
workflow {
    // Create input channel
    samples = Channel.fromPath("data/*.fastq")
        .map { file -> [file.baseName, file] }
    
    // Reference (same for all)
    reference = Channel.value("reference/hg38.fa")
    
    // Run QC
    qc = FASTQC(samples)
    
    // Align (uses reference)
    aligned = ALIGN(samples, reference)
    
    // View final results
    aligned.view()
}
```

</details>

---

## Exercise 4: Multi-Process Workflow (10 minutes)

### Part A: Understanding the Flow

**Given these processes:**

```groovy
process ALIGN {
    input: tuple val(id), path(fastq)
    output: tuple val(id), path("${id}.bam")
    script: """
        bwa mem reference.fa ${fastq} > ${id}.bam
    """
}

process CALL_VARIANTS {
    input: tuple val(id), path(bam)
    output: tuple val(id), path("${id}.vcf")
    script: """
        bcftools call ${bam} > ${id}.vcf
    """
}
```

**And this workflow:**

```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
        .map { f -> [f.baseName, f] }
    
    aligned = ALIGN(samples)
    variants = CALL_VARIANTS(aligned)
    variants.view()
}
```

**Answer:**
1. What channel does ALIGN output?
2. What does CALL_VARIANTS receive as input?
3. How many times does CALL_VARIANTS run?

<details>
<summary>✓ Solution A</summary>

**1. What channel does ALIGN output?**
```
Tuples like: [sample_id, /path/to/sample_id.bam]
Example: [tumor_001, /work/tumor_001.bam]
```

**2. What does CALL_VARIANTS receive as input?**
```
Each tuple from ALIGN's output:
- First call: [tumor_001, /work/tumor_001.bam]
- Second call: [tumor_002, /work/tumor_002.bam]
- Etc.
```

**3. How many times does CALL_VARIANTS run?**
```
Same number of times as input samples
If 5 FASTQ files → ALIGN runs 5 times → CALL_VARIANTS runs 5 times
All in parallel where possible!
```

**Data flow:**
```
samples channel (3 items)
    ↓
ALIGN (runs 3 times in parallel)
    ↓
aligned channel (3 items - BAM files)
    ↓
CALL_VARIANTS (runs 3 times in parallel)
    ↓
variants channel (3 items - VCF files)
    ↓
View results
```

</details>

---

### Part B: Add Filtering

**Modify the workflow to only process files larger than 1 GB:**

<details>
<summary>✓ Solution B</summary>

```groovy
workflow {
    // Create input and filter by size
    samples = Channel.fromPath("*.fastq")
        .filter { file -> file.size() > 1000000000 }  // > 1 GB
        .map { f -> [f.baseName, f] }
    
    aligned = ALIGN(samples)
    variants = CALL_VARIANTS(aligned)
    variants.view()
}
```

**Key addition:**
```groovy
.filter { file -> file.size() > 1000000000 }
```

This adds a filtering step after finding files but before processing.

**Timeline:**
```
Find all *.fastq files (maybe 10)
    ↓
Filter to only those > 1GB (maybe 7)
    ↓
Transform to tuples
    ↓
ALIGN (runs 7 times)
```

</details>

---

## Exercise 5: Real-World Pipeline (15 minutes)

### Challenge: Build a Complete QC Pipeline

**Requirements:**
- Input: Multiple FASTQ files in `data/`
- Step 1: Quality check with FastQC
- Step 2: Trim adapters with Trimmomatic
- Step 3: QC again on trimmed reads
- Output: Final QC reports

**Use these processes (already defined):**

```groovy
process FASTQC {
    input: tuple val(id), path(fastq)
    output: tuple val(id), path("${id}_raw_fastqc.html")
    script: """
        fastqc ${fastq}
        mv *_fastqc.html ${id}_raw_fastqc.html
    """
}

process TRIMMOMATIC {
    input: tuple val(id), path(fastq)
    output: tuple val(id), path("${id}.trimmed.fastq")
    script: """
        trimmomatic SE ${fastq} ${id}.trimmed.fastq ILLUMINACLIP:adapters.fa:2:30:10
    """
}

process FASTQC_POST {
    input: tuple val(id), path(fastq)
    output: tuple val(id), path("${id}_trimmed_fastqc.html")
    script: """
        fastqc ${fastq}
        mv *_fastqc.html ${id}_trimmed_fastqc.html
    """
}
```

**Your workflow:**

```groovy
workflow {
    ???
}
```

<details>
<summary>✓ Solution</summary>

```groovy
workflow {
    // Step 1: Create input channel
    samples = Channel.fromPath("data/*.fastq")
        .map { file -> [file.baseName, file] }
    
    // Step 2: QC on raw reads
    qc_raw = FASTQC(samples)
    
    // Step 3: Trim adapters
    trimmed = TRIMMOMATIC(samples)
    
    // Step 4: QC on trimmed reads
    qc_trimmed = FASTQC_POST(trimmed)
    
    // Step 5: View final results
    qc_trimmed.view()
}
```

**Data flow:**
```
samples channel
    ├─→ FASTQC ────→ qc_raw channel
    └─→ TRIMMOMATIC ────→ trimmed channel
                             ├─→ FASTQC_POST ────→ qc_trimmed channel
                             └─→ View results
```

**Key insight:**
- `samples` is used by both FASTQC and TRIMMOMATIC
- FASTQC and TRIMMOMATIC run in parallel on same input
- FASTQC_POST waits for trimmed output
- All 3 processes could run simultaneously on different samples!

**More production-ready version:**

```groovy
workflow {
    // Input configuration
    input_dir = "data/"
    
    // Step 1: Create and filter input
    samples = Channel.fromPath("${input_dir}/*.fastq")
        .filter { file -> file.size() > 100 }  // Ignore tiny files
        .map { file -> [file.baseName, file] }
    
    // Step 2: First QC
    qc_raw = FASTQC(samples)
    
    // Step 3: Trim
    trimmed = TRIMMOMATIC(samples)
    
    // Step 4: Second QC
    qc_trimmed = FASTQC_POST(trimmed)
    
    // Step 5: Publish results
    qc_trimmed.view()
    
    // (Optional) Collect all QC reports
    all_qc = qc_raw.mix(qc_trimmed).map { id, file -> file }.collect()
}
```

</details>

---

## Exercise 6: Debugging Workflows (10 minutes)

### Part A: Empty Workflow

**Your workflow produces no output:**

```groovy
workflow {
    samples = Channel.fromPath("data/*.fastq")
    sample_names = samples.map { f -> f.baseName }
    // No process call!
}
```

**What's wrong? How to fix?**

<details>
<summary>✓ Solution A</summary>

**Problem:**
- You created channels but never passed them to a process
- No process = no work = no output

**Fix:**
```groovy
workflow {
    samples = Channel.fromPath("data/*.fastq")
    sample_names = samples.map { f -> f.baseName }
    
    // Add process call
    MY_PROCESS(sample_names)  // Now it runs!
}
```

**Lesson:** Channels are just data in flight. You must pass them to processes to actually do work.

</details>

---

### Part B: Process Never Receives Data

**Your process never runs:**

```groovy
workflow {
    samples = Channel.fromPath("data/*.fastq")
    
    // Process defined elsewhere
    // But it's never called!
}
```

**What's wrong?**

<details>
<summary>✓ Solution B</summary>

**Problem:** 
Similar to Part A - channel is created but not passed to process.

**Fix:**
```groovy
workflow {
    samples = Channel.fromPath("data/*.fastq")
    
    // Pass channel to process
    results = MY_PROCESS(samples)
}
```

**Diagnostic tip:**
Run with `-resume` flag to see what's happening:
```bash
nextflow run workflow.nf -resume
```

</details>

---

## Exercise 7: Build Your First Real Workflow (20 minutes)

### Challenge: RNA-Seq Pipeline

**Build a workflow that:**
1. Finds paired-end FASTQ files
2. Runs FastQC
3. Aligns with STAR
4. Generates count matrix with featureCounts

**Processes provided:**

```groovy
process FASTQC {
    input: tuple val(id), path(reads)
    output: tuple val(id), path("${id}_fastqc.zip")
    script: """
        fastqc ${reads[0]} ${reads[1]}
    """
}

process ALIGN_STAR {
    cpus 8
    input:
        tuple val(id), path(reads)
        path star_index
    output:
        tuple val(id), path("Aligned.sortedByCoord.out.bam")
    script: """
        STAR --genomeDir ${star_index} \\
             --readFilesIn ${reads[0]} ${reads[1]} \\
             --runThreadN ${task.cpus} \\
             --outSAMtype BAM SortedByCoordinate
    """
}

process COUNT_FEATURES {
    input:
        tuple val(id), path(bam)
        path gtf_file
    output:
        tuple val(id), path("counts.txt")
    script: """
        featureCounts -a ${gtf_file} -o counts.txt ${bam}
    """
}
```

**Write your workflow:**

```groovy
workflow {
    ???
}
```

<details>
<summary>✓ Solution</summary>

```groovy
workflow {
    // Step 1: Find paired-end reads
    reads = Channel.fromFilePairs("data/*_{1,2}.fastq.gz")
    
    // Step 2: Static resources
    star_index = Channel.value("reference/STAR_index")
    gtf_file = Channel.value("reference/genes.gtf")
    
    // Step 3: QC (runs on both reads simultaneously)
    qc = FASTQC(reads)
    
    // Step 4: Alignment
    aligned = ALIGN_STAR(reads, star_index)
    
    // Step 5: Count features
    counts = COUNT_FEATURES(aligned, gtf_file)
    
    // Step 6: View and publish results
    counts.view()
}
```

**Explanation:**

1. **`Channel.fromFilePairs`** - Finds paired files
   - Matches sample1_1.fastq.gz and sample1_2.fastq.gz
   - Groups them as [sample1, [file1, file2]]

2. **`Channel.value`** - Static resources
   - STAR index used by all samples
   - GTF file used by all samples

3. **Parallel steps:**
   - FASTQC and ALIGN_STAR can run on different samples simultaneously
   - Once alignment done, COUNT_FEATURES can start

4. **Data preservation:**
   - Sample ID stays with data through entire pipeline
   - Easy to track what sample produced what counts

**Real-world notes:**
```groovy
workflow {
    // More robust version with error handling
    
    // Input configuration
    reads_pattern = "data/*_{1,2}.fastq.gz"
    star_index = "reference/STAR_index"
    gtf_file = "reference/genes.gtf"
    
    // Find reads (fail if not found)
    reads = Channel.fromFilePairs(reads_pattern)
        .ifEmpty { error("No paired reads found") }
    
    // QC
    qc = FASTQC(reads)
    
    // Align
    aligned = ALIGN_STAR(reads, Channel.value(star_index))
    
    // Count
    counts = COUNT_FEATURES(aligned, Channel.value(gtf_file))
    
    // Results
    counts.view()
    counts.map { id, file -> file }.collect().set { all_counts }
}
```

This version:
- Fails clearly if reads not found
- Better error messages
- Collects all counts if needed later

</details>

---

## Summary: What You've Practiced

✅ Creating channels (fromPath, from, value, fromFilePairs)  
✅ Transforming channels (map, filter)  
✅ Building simple workflows  
✅ Chaining processes  
✅ Understanding parallelization  
✅ Preserving sample metadata  
✅ Debugging workflow issues  
✅ Building real RNA-seq pipeline  

---

## 🎓 Self-Check

**Can you do these without looking at solutions?**

- [ ] Explain what a channel is
- [ ] Create a channel from files
- [ ] Transform a channel with .map()
- [ ] Filter a channel with .filter()
- [ ] Connect two processes in a workflow
- [ ] Preserve sample IDs through processes
- [ ] Debug a workflow that's not running

**If you can do most of these, you're ready for Day 5!**

---

## 🚀 What's Next (Day 5)

Tomorrow you'll learn about **subworkflows** and best practices for organizing larger pipelines.

You've learned:
- Processes ✅
- Workflows ✅
- Channels ✅

Next you'll learn:
- Organizing code
- Reusing workflows
- Best practices

Your first complete, production-ready pipeline is within reach!

---

*Day 4 of 28 - You've built your first complete Nextflow workflows! 🎉*
