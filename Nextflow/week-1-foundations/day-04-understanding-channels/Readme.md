# Day 4: Channels and Your First Complete Workflow

**Learning Time**: 30 minutes  
**Prerequisites**: Completed Days 1-3  
**Goal**: Build your first complete Nextflow workflow with multiple connected processes

---

## 📖 Introduction (2 minutes)

Welcome to Day 4! You've come a long way:
- **Day 1:** Learned what Nextflow is
- **Day 2:** Mastered Groovy syntax
- **Day 3:** Wrote your first process

Today is where everything comes together. You'll learn how to **connect processes into workflows** and watch data flow through your pipeline.

### The Missing Piece

You can write processes, but how do they talk to each other? How does data flow from one process to the next?

**Answer: Channels**

Channels are Nextflow's way of moving data between processes. Think of them as conveyor belts carrying data through your workflow.

---

## 🎯 Learning Objectives

By the end of Day 4, you should be able to:

1. **Explain channels** - What they are and why they matter
2. **Create channels** - From files, lists, values
3. **Transform channels** - Using map, filter, collect
4. **Connect processes** - Build workflows
5. **Write a complete workflow** - Multiple processes working together
6. **Understand data flow** - How data moves through your pipeline
7. **Debug workflows** - When things go wrong

---

## 📚 Key Concepts (20 minutes)

### 1. What is a Channel?

A **channel** is a stream of data that flows through your workflow.

Think of it like a conveyor belt:
```
Input files
    ↓
Channel: sample1.fastq → Process1 → sample1.bam → Channel → Process2 → sample1.vcf
Channel: sample2.fastq →   (parallel)  sample2.bam →   (parallel)   sample2.vcf
Channel: sample3.fastq →            sample3.bam →                sample3.vcf
```

**Key insight:** Channels enable parallelization. One channel carries multiple items, and Nextflow runs processes in parallel for each item.

### Types of Channels

#### Channel Type 1: Queue Channels
```groovy
channel = Channel.from(1, 2, 3, 4, 5)
// Emits: 1, 2, 3, 4, 5
// Each item goes to ONE downstream process
```

Used for:
- Passing data through workflows
- One-time use

#### Channel Type 2: Value Channels
```groovy
channel = Channel.value("reference.fa")
// Emits: "reference.fa"
// Available to ALL downstream processes
```

Used for:
- Reference genomes (needed by all)
- Parameters
- Configuration

---

### 2. Creating Channels: Common Patterns

#### Pattern 1: From Files
```groovy
// Single files
reads = Channel.fromPath("data/sample.fastq")

// Multiple files (most common!)
samples = Channel.fromPath("data/*.fastq")
```

When `*.fastq` matches 3 files:
- Channel emits 3 items
- Processes run 3 times (in parallel)

#### Pattern 2: From a List
```groovy
// Create channel from a list
samples = Channel.from("sample1", "sample2", "sample3")

// Or from actual list
sample_list = ["s1", "s2", "s3"]
samples = Channel.from(sample_list)
```

#### Pattern 3: From Pairs (Paired-End Reads)
```groovy
// Most common for sequencing data
reads = Channel.fromFilePairs("data/*_{1,2}.fastq.gz")
// Emits: ["sample1", [/path/R1.fastq, /path/R2.fastq]]
```

The `{1,2}` pattern matches:
- `sample1_1.fastq` and `sample1_2.fastq`
- Creates tuples keeping pairs together

#### Pattern 4: From CSV
```groovy
samples = Channel.fromPath("metadata.csv")
    .splitCsv(header:true)
    // Now each line is a map: [id: "s1", condition: "control", ...]
```

---

### 3. Transforming Channels: The Power of Groovy

Channels are streams, and you transform them with closures (remember Day 2?).

#### Transform 1: Map (Transform Each Item)
```groovy
// Input channel has file paths
files = Channel.fromPath("data/*.fastq")

// Transform to sample names
names = files.map { file -> file.baseName }

// Or more complex transformation
samples = files.map { file ->
    [file.baseName, file]  // Returns tuple
}
```

**Example output:**
```
Input:  /data/tumor_001.fastq
Output: tumor_001

Input:  [/data/tumor_001.fastq]
Output: [tumor_001, /data/tumor_001.fastq]
```

#### Transform 2: Filter (Keep Matching Items)
```groovy
// Input: files with various sizes
files = Channel.fromPath("data/*.fastq")

// Keep only files > 1GB
large_files = files.filter { file -> file.size() > 1000000000 }
```

#### Transform 3: Collect (Gather All Items)
```groovy
// Input: multiple items in channel
samples = Channel.from(1, 2, 3, 4, 5)

// Collect all into one list
all_samples = samples.collect()
// Output: [1, 2, 3, 4, 5] as a single item

// Practical example: collect BAM files
bam_files = Channel.fromPath("results/*.bam")
all_bams = bam_files.collect()
// Now pass this to a merge process that needs all files at once
```

**Critical difference:**
- Without `collect()`: Pass items one at a time (parallelized)
- With `collect()`: Wait for all items, pass together (single run)

#### Transform 4: Join (Combine Channels)
```groovy
// Two channels with matching keys
reads = Channel.fromPath("data/*.fastq").map { f -> [f.baseName, f] }
metadata = Channel.from([
    ["s1", "control"],
    ["s2", "treatment"]
])

// Join them
combined = reads.join(metadata)
// Output: [s1, /path/s1.fastq, control]
```

#### Transform 5: Branch (Split Channel)
```groovy
// Input channel with tuples [id, file, condition]
samples = Channel.from(
    ["s1", "file1", "control"],
    ["s2", "file2", "treatment"],
    ["s3", "file3", "control"]
)

// Split into two channels
branched = samples.branch {
    control: it[2] == "control"
    treatment: it[2] == "treatment"
}

// Now use separately
branched.control.view()     // Shows control samples
branched.treatment.view()   // Shows treatment samples
```

---

### 4. Building a Workflow: The Complete Picture

A **workflow** is where you:
1. Create channels
2. Transform channels (optionally)
3. Connect processes together

```groovy
workflow {
    // Step 1: Create input channel
    samples = Channel.fromPath("data/*.fastq")
    
    // Step 2: Transform (optional)
    // samples = samples.filter { it.size() > 100 }
    
    // Step 3: Run first process
    qc_results = QUALITY_CHECK(samples)
    
    // Step 4: Chain to next process
    aligned = ALIGNMENT(samples)
    
    // Step 5: Chain to third process
    variants = VARIANT_CALLING(aligned)
    
    // Step 6: View final results
    variants.view()
}
```

**The flow:**
```
samples channel
    ↓
QUALITY_CHECK process (runs for each sample in parallel)
    ↓
qc_results channel
    ↓
ALIGNMENT process (runs for each sample in parallel)
    ↓
aligned channel
    ↓
VARIANT_CALLING process (runs for each sample in parallel)
    ↓
variants channel
    ↓
View results
```

---

### 5. Real Example: Complete Quality Control Workflow

Let's build a real workflow with 3 processes:

```groovy
// Define processes (from previous days)
process FASTQC {
    cpus 2
    input:
        path fastq
    output:
        path "*_fastqc.html"
    script:
        """
        fastqc ${fastq}
        """
}

process ADAPTER_TRIM {
    cpus 4
    input:
        path fastq
    output:
        path "*.trimmed.fastq"
    script:
        """
        trim_galore ${fastq}
        """
}

process SECOND_QC {
    cpus 2
    input:
        path fastq
    output:
        path "*_fastqc.html"
    script:
        """
        fastqc ${fastq}
        """
}

// Now define the workflow
workflow {
    // Create input channel
    raw_reads = Channel.fromPath("data/*.fastq")
    
    // Run initial QC
    qc1 = FASTQC(raw_reads)
    
    // Trim adapters
    trimmed = ADAPTER_TRIM(raw_reads)
    
    // Run QC on trimmed reads
    qc2 = SECOND_QC(trimmed)
    
    // View final reports
    qc2.view()
}
```

**What happens:**
1. Find all FASTQ files in `data/`
2. Run FastQC on all of them (in parallel)
3. Trim adapters on all of them (in parallel)
4. Run QC on trimmed reads
5. All processes can run simultaneously!

---

### 6. Practical: Multi-Tool Pipeline

A more realistic workflow with real tools:

```groovy
process QUALITY_CHECK {
    input:
        tuple val(sample_id), path(fastq)
    output:
        tuple val(sample_id), path("${sample_id}_fastqc.html")
    script:
        """
        fastqc ${fastq}
        mv *_fastqc.html ${sample_id}_fastqc.html
        """
}

process ALIGNMENT {
    cpus 8
    input:
        tuple val(sample_id), path(fastq)
        path reference
    output:
        tuple val(sample_id), path("${sample_id}.bam")
    script:
        """
        bwa mem -t ${task.cpus} ${reference} ${fastq} | samtools view -b > ${sample_id}.bam
        samtools sort -o ${sample_id}.bam ${sample_id}.bam
        """
}

process VARIANT_CALLING {
    input:
        tuple val(sample_id), path(bam)
        path reference
    output:
        tuple val(sample_id), path("${sample_id}.vcf")
    script:
        """
        bcftools mpileup -f ${reference} ${bam} | bcftools call -mv -o ${sample_id}.vcf
        """
}

// The workflow
workflow {
    // Input: multiple FASTQ files
    samples = Channel.fromPath("data/*.fastq")
        .map { file -> [file.baseName, file] }  // Create tuples with sample names
    
    // Reference genome (available to all processes)
    reference = Channel.value("reference/hg38.fa")
    
    // Step 1: Quality check (input: sample tuples)
    qc = QUALITY_CHECK(samples)
    
    // Step 2: Alignment (input: sample tuples + reference)
    aligned = ALIGNMENT(samples, reference)
    
    // Step 3: Variant calling (input: aligned BAMs + reference)
    variants = VARIANT_CALLING(aligned, reference)
    
    // View results
    variants.view()
}
```

**Key patterns here:**
- Input files are turned into tuples: `[sample_id, file]`
- Sample ID travels with the file through all processes
- Reference is a value channel (available to all)
- Each process returns output in same tuple format

---

### 7. Understanding: How Nextflow Schedules Work

When you have this workflow:

```groovy
workflow {
    samples = Channel.fromPath("*.fastq")  // 3 files
    
    aligned = ALIGNMENT(samples)     // Align each
    variants = VARIANT_CALLING(aligned)  // Call variants on each
}
```

**Timeline (assuming 4 cores available):**
```
Time 0s:
  ALIGNMENT_s1 (core 1) ████████████
  ALIGNMENT_s2 (core 2) ████████████
  ALIGNMENT_s3 (core 3) ████████████

Time 10s (after first alignment completes):
  ALIGNMENT_s4 (core 1) ████████████
  VARIANT_s1 (core 2)   █████████
  VARIANT_s2 (core 3)   █████████

Time 20s (after first variant call):
  VARIANT_s3 (core 1)   █████████
  VARIANT_s4 (core 2)   █████████
```

**Key insight:** Nextflow schedules smartly. As soon as one process produces output, the next process can start if data is available. This maximizes parallelization.

---

### 8. Common Channel Patterns

#### Pattern 1: Single Input File
```groovy
workflow {
    reference = Channel.value("reference.fa")
    PROCESS(reference)
}
```

#### Pattern 2: Multiple Input Files
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    PROCESS(samples)  // Runs for each file
}
```

#### Pattern 3: Reference + Input
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    reference = Channel.value("ref.fa")
    PROCESS(samples, reference)  // Ref used by all
}
```

#### Pattern 4: Paired Data
```groovy
workflow {
    reads = Channel.fromFilePairs("*_{1,2}.fastq.gz")
    PROCESS(reads)  // Each pair goes as unit
}
```

#### Pattern 5: Multiple Outputs to Downstream
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    
    aligned = ALIGN(samples)
    qc = QUALITY_CHECK(samples)  // Both use same input
    
    // aligned and qc are two separate channels
}
```

---

### 9. Debugging Workflows: Common Issues

#### Issue 1: Channel Doesn't Emit Any Items
```groovy
samples = Channel.fromPath("data/*.fastq")  // No .fastq files found!
```
**Fix:** Check file path and pattern match

#### Issue 2: Process Never Runs
```groovy
workflow {
    samples = Channel.from("s1", "s2")
    // Created channel but never passed to process!
    // Need: PROCESS(samples)
}
```

#### Issue 3: Wrong Number of Outputs
```groovy
process MY_PROCESS {
    output:
        path "*.txt"
    script:
        """
        echo "test" > output1.txt
        # Missing: echo "test2" > output2.txt
        """
}
```
**Fix:** Make sure script creates files matching output pattern

---

### 10. The Complete Workflow Template

```groovy
// 1. Define your processes (from Day 3)
process STEP1 {
    input:
        path input_file
    output:
        path "output.txt"
    script:
        """
        command ${input_file} > output.txt
        """
}

process STEP2 {
    input:
        path input_file
    output:
        path "result.txt"
    script:
        """
        command2 ${input_file} > result.txt
        """
}

// 2. Define your workflow
workflow {
    // Create input channels
    input_data = Channel.fromPath("data/*.fastq")
    
    // (Optional) Transform channels
    // input_data = input_data.filter { it.size() > 100 }
    
    // Chain processes
    step1_out = STEP1(input_data)
    final_out = STEP2(step1_out)
    
    // (Optional) View results
    final_out.view()
}

// 3. Run it
// $ nextflow run workflow.nf
```

---

## 🔗 Connecting to Day 3

### Recap: Day 3 Gave You Processes

```groovy
process FASTQC {
    input: path fastq
    output: path "*.html"
    script: "fastqc ${fastq}"
}
```

### Today: Day 4 Gives You Workflows

```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    FASTQC(samples)
}
```

**Together:** You can now run bioinformatics pipelines!

---

## ✅ Completion Checklist

- [ ] I understand what channels are
- [ ] I can create channels from files
- [ ] I can use .map() to transform data
- [ ] I can use .filter() to select items
- [ ] I know when to use .collect()
- [ ] I can connect processes in a workflow
- [ ] I understand parallelization
- [ ] I can debug workflow issues
- [ ] I feel ready to build workflows

---

## 🔑 Key Takeaways

### What is a Channel?
A stream of data flowing through your workflow. Multiple items = automatic parallelization.

### How Do I Create Channels?
```groovy
Channel.fromPath("*.fastq")      // From files
Channel.from(1, 2, 3)            // From list
Channel.value("reference.fa")    // Single value
Channel.fromFilePairs("*_{1,2}") // Paired files
```

### How Do I Transform Channels?
```groovy
.map { item -> transform(item) }           // Transform
.filter { item -> condition(item) }        // Filter
.collect()                                  // Gather all
.join(other_channel)                       // Combine
```

### How Do I Build a Workflow?
1. Create channels
2. Transform (optional)
3. Pass to processes
4. Chain processes together

---

## 🚀 Ready for Exercises?

You now understand:
- Channels ✅
- Workflows ✅
- Data flow ✅
- Process connections ✅

Time to build your first complete workflow!

---

*This is Day 4 of 28. You're building your first real bioinformatics pipelines!*
