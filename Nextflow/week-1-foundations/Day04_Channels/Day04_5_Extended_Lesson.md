# Day 4: Understanding Channels

**Learning Time**: 30 minutes  
**Prerequisites**: Days 1-3 completed, understanding of processes  
**Goal**: Master channels as the data highways of Nextflow workflows

---

## 📖 Introduction (3 minutes)

Welcome to Day 4! Today you'll learn about **channels**—the most conceptually unique part of Nextflow. This is where Nextflow truly differs from traditional programming, and understanding channels is the key to unlocking Nextflow's power.

**Why channels matter**: They're the mechanism that enables automatic parallelization. When you understand channels, you understand how Nextflow can process 1,000 samples as easily as it processes 1.

### What You'll Learn Today

- What channels actually are (not just "data pipes")
- The crucial difference between queue and value channels
- How to create channels from files, lists, and values
- How channels connect to processes
- Why channels enable automatic parallelization
- Common patterns for working with channel data

### The Big Insight

In Python, you load data, then loop through it:
```python
files = glob.glob("*.fastq")  # Load all at once
for file in files:             # Explicit loop
    process(file)              # Sequential
```

In Nextflow, data **flows** through channels and processes:
```groovy
Channel.fromPath("*.fastq")    // Create a stream
    | processFile              // Data flows through automatically
```

Channels are **not** variables holding data. They're **streams** that emit data to processes.

---

## 🎯 Key Concepts with Examples (12 minutes)

### 1. What Is a Channel? (2 minutes)

**Conceptual Definition**: A channel is a **data stream** that emits values, one at a time, to processes that consume them.

Think of a channel as a **conveyor belt in a factory**:
- Items (data) are placed on the belt
- Items move along the belt one by one
- Workers (processes) take items as they arrive
- Multiple workers can work on different items simultaneously

**Key Insight**: When a process receives input from a channel with multiple items, Nextflow automatically runs the process multiple times—once for each item, in parallel.

```groovy
// This channel emits 3 items
Channel.from("sample1", "sample2", "sample3")

// When connected to a process, that process runs 3 times
// (potentially in parallel!)
```

### 2. Channel Types: Queue vs Value (3 minutes)

Nextflow has two fundamental channel types:

#### **Queue Channels** (Default)

- Can be consumed **only once**
- Emit multiple values
- Used for data that flows through the pipeline
- When all values are emitted, the channel closes

```groovy
// Queue channel with 3 values
samples = Channel.from("A", "B", "C")

process1(samples)  // Consumes: A, B, C
process2(samples)  // ❌ ERROR! Already consumed!
```

**Think**: Like items on a conveyor belt—once taken, they're gone.

#### **Value Channels**

- Can be consumed **multiple times**
- Emit a **single value** (or tuple)
- Used for reference data that many processes need
- Never closes

```groovy
// Value channel with single value
reference = Channel.value("reference_genome.fa")

process1(reference)  // Uses reference
process2(reference)  // ✅ Can reuse!
process3(reference)  // ✅ Can reuse!
```

**Think**: Like a reference manual that everyone can read without consuming it.

#### **When to Use Each**

| Scenario | Channel Type | Example |
|----------|--------------|---------|
| Multiple samples to process | Queue | `Channel.fromPath("*.fastq")` |
| Reference genome for all samples | Value | `Channel.value(params.genome)` |
| List of parameters | Queue | `Channel.from(1..10)` |
| Single configuration value | Value | `Channel.value(params.threads)` |

**Critical Rule**: If you need to use the same data in multiple processes, use a **value channel** or explicitly split a queue channel.

### 3. Creating Channels (3 minutes)

Nextflow provides several channel factories:

#### **fromPath** - Files and Directories
```groovy
// Single pattern
fastq_files = Channel.fromPath("data/*.fastq")
// Emits: data/sample1.fastq, data/sample2.fastq, ...

// Multiple patterns
all_files = Channel.fromPath("data/*.{fastq,fq}")
// Emits: files matching either extension

// Recursive search
deep_files = Channel.fromPath("data/**/*.fastq")
// Emits: all .fastq files in data/ and subdirectories

// With options
files = Channel.fromPath("data/*.fastq", checkIfExists: true)
// Throws error if no files match
```

#### **fromFilePairs** - Paired Files
```groovy
// For paired-end reads: sample_R1.fastq, sample_R2.fastq
paired_reads = Channel.fromFilePairs("data/*_{R1,R2}.fastq")
// Emits: [sample, [sample_R1.fastq, sample_R2.fastq]]

// More complex pattern
paired_reads = Channel.fromFilePairs("data/*_{1,2}.fastq.gz")
// Emits tuples: [sample_id, [read1, read2]]
```

#### **from** / **of** - Values and Lists
```groovy
// From individual values
numbers = Channel.from(1, 2, 3, 4, 5)
// Emits: 1, 2, 3, 4, 5 (one at a time)

// From a list
samples = Channel.from(["A", "B", "C"])
// Emits: A, B, C

// From a range
range = Channel.from(1..100)
// Emits: 1, 2, 3, ..., 100

// Modern syntax (Nextflow 20.07+)
numbers = Channel.of(1, 2, 3, 4, 5)  // Same as .from()
```

#### **value** - Single Reusable Value
```groovy
// Reference genome (used by all samples)
genome = Channel.value("/data/reference/hg38.fa")

// Number of threads (configuration value)
threads = Channel.value(params.threads)

// Any single value you need to reuse
constant = Channel.value("this value can be used multiple times")
```

#### **fromList** - Create from Groovy List
```groovy
// Create from existing Groovy list
sample_list = ["sample1", "sample2", "sample3"]
channel = Channel.fromList(sample_list)
// Emits: sample1, sample2, sample3
```

### 4. Channels with Processes (2 minutes)

This is where it all comes together! Let's see how channels feed processes:

```groovy
// Define a process (from Day 3)
process countReads {
    input:
    path fastq
    
    output:
    path "count.txt"
    
    script:
    """
    wc -l ${fastq} > count.txt
    """
}

// Create a channel
fastq_ch = Channel.fromPath("data/*.fastq")

// Connect channel to process
workflow {
    counts = countReads(fastq_ch)
}
```

**What happens**:
1. Channel emits each .fastq file
2. For each file, Nextflow runs `countReads` once
3. If there are 100 files, `countReads` runs 100 times
4. All runs happen in parallel (respecting resource limits)
5. Each run produces its own `count.txt`
6. All outputs are collected into a new channel: `counts`

**The Magic**: You write the process once, Nextflow handles the parallelization!

#### **Process with Multiple Inputs**

```groovy
process align {
    input:
    tuple val(sample_id), path(reads)
    path reference  // Value channel - same for all
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    script:
    """
    bwa mem ${reference} ${reads} > ${sample_id}.bam
    """
}

workflow {
    // Queue channel: different for each sample
    samples = Channel.fromFilePairs("data/*_{1,2}.fastq")
    
    // Value channel: same reference for all
    ref = Channel.value("genome.fa")
    
    // Each sample uses the same reference
    aligned = align(samples, ref)
}
```

**Understanding**:
- `samples` is a queue channel → different data for each process run
- `ref` is a value channel → same data for all process runs
- Nextflow runs `align` once per sample, all using the same reference

### 5. Automatic Parallelization Explained (2 minutes)

This is the killer feature. Let's see it in action:

**Scenario**: 100 FASTQ files need quality control

```groovy
process fastqc {
    input:
    path fastq
    
    output:
    path "*.html"
    
    script:
    """
    fastqc ${fastq}
    """
}

workflow {
    fastq_files = Channel.fromPath("samples/*.fastq")
    // Suppose this emits 100 files
    
    qc_results = fastqc(fastq_files)
    // Nextflow automatically:
    // 1. Creates 100 process instances
    // 2. Distributes them across available CPUs
    // 3. Runs as many in parallel as resources allow
    // 4. Queues the rest until resources are available
}
```

**Time Comparison**:

| Approach | Time for 100 Samples (5 min each) |
|----------|-----------------------------------|
| Python sequential loop | 500 minutes (8.3 hours) |
| Nextflow with 8 cores | ~65 minutes |
| Nextflow with 32 cores | ~20 minutes |
| Nextflow on cluster (unlimited) | ~5 minutes |

**You wrote the same code**, Nextflow scaled it automatically!

---

## 🔗 Python Connection: Iterables vs Channels (3 minutes)

### Python: Eager Loading and Explicit Loops

```python
import glob

# Load all files into memory at once
files = glob.glob("data/*.fastq")  # List with all filenames

# Explicit sequential loop
for file in files:
    result = process_file(file)
    results.append(result)

# For parallel processing, you need:
from concurrent.futures import ProcessPoolExecutor

with ProcessPoolExecutor(max_workers=8) as executor:
    # Manually manage parallelization
    results = list(executor.map(process_file, files))
```

**Characteristics**:
- ✅ Familiar, explicit control
- ❌ All data loaded at once (memory issue for large datasets)
- ❌ Sequential by default
- ❌ Manual parallelization code needed
- ❌ Hard to scale beyond one machine

### Nextflow: Lazy Streaming and Implicit Parallelization

```groovy
// Create a stream of files (not all loaded at once)
files = Channel.fromPath("data/*.fastq")

// Process automatically parallelizes
workflow {
    results = process_file(files)
    // Nextflow handles:
    // - How many to run in parallel
    // - When to schedule each
    // - Resource allocation
    // - Distribution across cluster nodes
}
```

**Characteristics**:
- ✅ Streaming (low memory footprint)
- ✅ Parallel by default
- ✅ Automatic parallelization
- ✅ Scales from laptop to cluster with same code
- ⚠️ Less explicit control (different paradigm)

### Key Conceptual Differences

| Concept | Python List/Iterator | Nextflow Channel |
|---------|---------------------|------------------|
| **Loading** | Eager (load all) | Lazy (stream) |
| **Consumption** | Multiple times | Once (queue) or unlimited (value) |
| **Parallelization** | Manual | Automatic |
| **Memory** | All items in memory | One at a time |
| **Control flow** | Explicit loops | Implicit (data-driven) |
| **Scaling** | Single machine | Multi-machine ready |

### Mental Model Shift

**Python thinking**: "Load the data, then loop through it"
```python
data = load_all()        # Load
for item in data:        # Loop
    process(item)        # Process
```

**Nextflow thinking**: "Create a stream, let it flow through"
```groovy
data = Channel.create()  // Stream
data | process          // Flow through (automatic loop)
```

### The List Analogy (But It's Wrong!)

People often think: "Channel is like a Python list"

**Why this is misleading**:
```python
# Python list
my_list = [1, 2, 3]
function1(my_list)  # List still exists
function2(my_list)  # Can reuse
```

```groovy
// Nextflow queue channel
my_channel = Channel.from(1, 2, 3)
process1(my_channel)  // Channel consumed
process2(my_channel)  // ❌ ERROR! Empty!
```

**Better analogy**: Channel is like a **Python generator**
```python
# Python generator (similar to queue channel)
def my_generator():
    yield 1
    yield 2
    yield 3

gen = my_generator()
list(gen)  # [1, 2, 3]
list(gen)  # [] - exhausted!
```

But even this is imperfect because Nextflow channels enable **automatic parallelization** which Python generators don't.

---

## 💻 Hands-On Exercises (10 minutes)

### Exercise 1: Creating Different Channel Types (3 minutes)

**Task**: Create channels for different data scenarios.

**Scenario A**: You have FASTQ files named `sample1.fastq`, `sample2.fastq`, ... in a directory `data/`.

Create a channel that emits each file.

<details>
<summary>Click for solution</summary>

```groovy
// Solution A
fastq_ch = Channel.fromPath("data/*.fastq")

// Alternative with checking
fastq_ch = Channel.fromPath("data/*.fastq", checkIfExists: true)

// If in subdirectories too
fastq_ch = Channel.fromPath("data/**/*.fastq")
```

**What it emits**: Each matching file path as it's discovered
```
data/sample1.fastq
data/sample2.fastq
data/sample3.fastq
...
```
</details>

**Scenario B**: You have a reference genome `hg38.fa` that every sample needs to use.

Create a channel that can be used by multiple processes.

<details>
<summary>Click for solution</summary>

```groovy
// Solution B
reference = Channel.value("hg38.fa")

// Or with full path
reference = Channel.value("/data/references/hg38.fa")

// Or from a parameter
reference = Channel.value(params.genome)

// Usage in workflow
workflow {
    samples = Channel.fromPath("*.fastq")
    
    process1(samples, reference)  // Uses reference
    process2(samples, reference)  // Reuses reference ✓
    process3(samples, reference)  // Reuses reference ✓
}
```

**Why value channel**: The same reference is used repeatedly, so it must be a value channel (reusable).
</details>

**Scenario C**: You have paired-end reads: `sample1_R1.fastq`, `sample1_R2.fastq`, etc.

Create a channel that pairs them correctly.

<details>
<summary>Click for solution</summary>

```groovy
// Solution C
paired_ch = Channel.fromFilePairs("data/*_{R1,R2}.fastq")

// More flexible pattern
paired_ch = Channel.fromFilePairs("data/*_{1,2}.fastq.gz")

// With custom grouping
paired_ch = Channel.fromFilePairs("data/*_R{1,2}.fastq") {
    file -> file.name.replaceAll(/_R[12].fastq$/, '')
}
```

**What it emits**: Tuples of [sample_id, [read1, read2]]
```
[sample1, [sample1_R1.fastq, sample1_R2.fastq]]
[sample2, [sample2_R1.fastq, sample2_R2.fastq]]
[sample3, [sample3_R1.fastq, sample3_R2.fastq]]
...
```

**In a process**:
```groovy
process alignPaired {
    input:
    tuple val(sample_id), path(reads)
    
    script:
    """
    bwa mem ref.fa ${reads[0]} ${reads[1]} > ${sample_id}.sam
    # Or: bwa mem ref.fa ${reads} (works too!)
    """
}
```
</details>

**Scenario D**: You need to test different quality thresholds: 20, 25, 30, 35, 40.

Create a channel with these values.

<details>
<summary>Click for solution</summary>

```groovy
// Solution D - Multiple ways:

// Method 1: Explicit values
thresholds = Channel.from(20, 25, 30, 35, 40)

// Method 2: From a list
threshold_list = [20, 25, 30, 35, 40]
thresholds = Channel.fromList(threshold_list)

// Method 3: Using a range with step
thresholds = Channel.from(20, 25, 30, 35, 40)

// Method 4: If they were sequential
thresholds = Channel.from(20..40)  // 20, 21, 22, ..., 40

// Method 5: Modern syntax
thresholds = Channel.of(20, 25, 30, 35, 40)
```

**Usage**:
```groovy
process filterByQuality {
    input:
    path vcf
    each threshold  // Combines vcf × each threshold
    
    output:
    path "filtered_q${threshold}.vcf"
    
    script:
    """
    bcftools filter -i "QUAL>${threshold}" ${vcf} > filtered_q${threshold}.vcf
    """
}

workflow {
    vcf = Channel.value("variants.vcf")
    thresholds = Channel.from(20, 25, 30, 35, 40)
    
    filterByQuality(vcf, thresholds)
    // Runs 5 times, once per threshold
}
```
</details>

### Exercise 2: Understanding Queue vs Value Channels (3 minutes)

**Task**: Predict what happens in each scenario.

**Scenario A**:
```groovy
process printSample {
    input:
    val sample
    
    script:
    """
    echo "Processing ${sample}"
    """
}

workflow {
    samples = Channel.from("A", "B", "C")
    
    printSample(samples)
    printSample(samples)  // What happens here?
}
```

**Question**: Will this work? What will be printed?

<details>
<summary>Click for answer</summary>

**Answer**: ❌ This will **fail** with an error!

**Why**: `samples` is a **queue channel**. After the first `printSample(samples)`, the channel is consumed (empty). The second `printSample(samples)` tries to read from an empty channel.

**Error message**: Something like "Channel has been used as input by more than one process or operator"

**How to fix**:

**Option 1**: Use a value channel (if only one sample)
```groovy
sample = Channel.value("A")
printSample(sample)
printSample(sample)  // ✓ Works
```

**Option 2**: Split the queue channel
```groovy
samples = Channel.from("A", "B", "C")
samples.into { samples1; samples2 }  // Deprecated in newer Nextflow

// Modern way:
samples = Channel.from("A", "B", "C")
samples.multiMap { sample ->
    first: sample
    second: sample
}.set { result }

printSample(result.first)
printSample(result.second)
```

**Option 3**: Call both processes with the same channel reference (only works if they're in sequence)
```groovy
samples = Channel.from("A", "B", "C")
result1 = printSample(samples)
// Use result1 as input to next process, not original samples
```
</details>

**Scenario B**:
```groovy
workflow {
    genome = Channel.value("hg38.fa")
    
    process1(genome)
    process2(genome)
    process3(genome)
}
```

**Question**: Will this work? How many times will each process run?

<details>
<summary>Click for answer</summary>

**Answer**: ✅ This **works perfectly**!

**Why**: `genome` is a **value channel**. It can be consumed unlimited times.

**How many times each runs**: Each process runs **once**, receiving "hg38.fa" as input.

**Output**:
- process1 runs 1 time with "hg38.fa"
- process2 runs 1 time with "hg38.fa"
- process3 runs 1 time with "hg38.fa"

This is the correct pattern for reference data that multiple processes need.
</details>

### Exercise 3: Connecting Channels to Processes (4 minutes)

**Task**: Complete this workflow.

You have:
- FASTQ files in `data/` directory
- A reference genome file at `refs/genome.fa`
- A process `alignReads` that takes a FASTQ file and reference genome

**Complete this workflow**:
```groovy
process alignReads {
    input:
    path fastq
    path reference
    
    output:
    path "aligned.bam"
    
    script:
    """
    bwa mem ${reference} ${fastq} | samtools view -b > aligned.bam
    """
}

workflow {
    // TODO: Create channel for FASTQ files
    
    // TODO: Create channel for reference genome
    
    // TODO: Connect to alignReads process
}
```

<details>
<summary>Click for solution</summary>

```groovy
process alignReads {
    input:
    path fastq
    path reference
    
    output:
    path "aligned.bam"
    
    script:
    """
    bwa mem ${reference} ${fastq} | samtools view -b > aligned.bam
    """
}

workflow {
    // Create queue channel for FASTQ files (multiple samples)
    fastq_ch = Channel.fromPath("data/*.fastq")
    
    // Create value channel for reference (single, reused)
    reference_ch = Channel.value("refs/genome.fa")
    
    // Connect to process
    alignReads(fastq_ch, reference_ch)
    
    // Or store results:
    bam_files = alignReads(fastq_ch, reference_ch)
}
```

**What happens**:
1. `fastq_ch` emits each .fastq file
2. For each FASTQ file:
   - `alignReads` runs once
   - It receives the FASTQ file and the reference genome
   - It produces a BAM file
3. If there are 10 FASTQ files, `alignReads` runs 10 times in parallel
4. All 10 runs use the same reference genome (from value channel)

**Why these channel types**:
- **Queue channel for FASTQ**: Different file for each run
- **Value channel for reference**: Same file for all runs

**Better version with sample tracking**:
```groovy
process alignReads {
    input:
    tuple val(sample_id), path(fastq)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    script:
    """
    bwa mem ${reference} ${fastq} | samtools view -b > ${sample_id}.bam
    """
}

workflow {
    // Create channel with sample IDs
    fastq_ch = Channel.fromPath("data/*.fastq")
        .map { file -> [file.baseName, file] }
    
    reference_ch = Channel.value("refs/genome.fa")
    
    alignReads(fastq_ch, reference_ch)
}
```

This keeps sample IDs with their BAM files for downstream processing!
</details>

---

## 🤔 Reflection Activity (4 minutes)

### Question 1: The Parallelization Insight

Consider this Python code:
```python
samples = ["s1", "s2", "s3", "s4", "s5"]
results = []

for sample in samples:
    result = expensive_analysis(sample)  # Takes 1 hour
    results.append(result)

# Total time: 5 hours
```

**Questions**:
1. How would you make this parallel in Python?
2. How does Nextflow handle this automatically?
3. What if you had 1,000 samples instead of 5?

<details>
<summary>Click for thoughts</summary>

**1. Python parallelization**:
```python
from concurrent.futures import ProcessPoolExecutor

samples = ["s1", "s2", "s3", "s4", "s5"]

with ProcessPoolExecutor(max_workers=5) as executor:
    results = list(executor.map(expensive_analysis, samples))

# Time: ~1 hour (if 5+ cores available)
```

You need to:
- Import parallelization libraries
- Manage worker pools
- Handle result collection
- Deal with exceptions in parallel contexts
- Worry about resource limits

**2. Nextflow automatic handling**:
```groovy
process expensiveAnalysis {
    input:
    val sample
    
    output:
    path "result.txt"
    
    script:
    """
    expensive_analysis ${sample} > result.txt
    """
}

workflow {
    samples = Channel.from("s1", "s2", "s3", "s4", "s5")
    results = expensiveAnalysis(samples)
}

// Time: ~1 hour (automatically parallel)
```

Nextflow:
- Automatically parallelizes across available resources
- Handles resource management
- Manages work directories
- Provides automatic retry on failure
- Scales from laptop to cluster with no code changes

**3. With 1,000 samples**:

**Python**: You'd need to carefully tune worker pool size, manage memory, potentially batch the work.

**Nextflow**: **Same code!** It automatically:
- Schedules tasks based on available resources
- Queues tasks that can't run immediately
- Distributes across cluster nodes if available
- Manages thousands of concurrent tasks effortlessly

The code doesn't change from 5 to 1,000 to 10,000 samples!
</details>

### Question 2: Queue vs Value Channel Decision

For each scenario, decide: **Queue channel** or **Value channel**?

**A.** FASTQ files to analyze (100 files)  
**B.** Database file for BLAST searches (1 file, used by all)  
**C.** Sample metadata CSV file (1 file, used by all)  
**D.** List of chromosomes to process (chr1, chr2, ..., chr22)  
**E.** Output directory path (same for all processes)  

<details>
<summary>Click for answers</summary>

**A. Queue channel** ✅
```groovy
fastq_ch = Channel.fromPath("data/*.fastq")
```
- 100 different files
- Each processed separately
- Queue channel emits each one

**B. Value channel** ✅
```groovy
database = Channel.value("blast_db/nr.fa")
```
- Single file
- Used by many samples
- Needs to be reusable

**C. Value channel** ✅
```groovy
metadata = Channel.value("sample_metadata.csv")
```
- Single file
- Shared reference
- Multiple processes might need it

**D. Queue channel** ✅
```groovy
chromosomes = Channel.from("chr1", "chr2", ..., "chr22")
// Or
chromosomes = Channel.from(1..22).map { "chr${it}" }
```
- Different values
- Process each separately
- Queue channel emits each chromosome

**E. Value channel** ✅
```groovy
outdir = Channel.value(params.output_dir)
// Or just use params.output_dir directly in processes
```
- Single value
- Used consistently
- Configuration parameter

**The Pattern**:
- **Multiple different things** → Queue channel
- **One thing used repeatedly** → Value channel
</details>

### Question 3: Visualizing Data Flow

Draw (mentally or on paper) the data flow for this workflow:

```groovy
process trim {
    input: path(fastq)
    output: path("trimmed.fq")
    script: "trim_tool ${fastq} > trimmed.fq"
}

process align {
    input:
    path(trimmed)
    path(reference)
    output: path("aligned.bam")
    script: "bwa mem ${reference} ${trimmed} | samtools view -b > aligned.bam"
}

workflow {
    reads = Channel.fromPath("data/*.fastq")  // 3 files
    ref = Channel.value("genome.fa")
    
    trimmed = trim(reads)
    aligned = align(trimmed, ref)
}
```

**Questions**:
1. How many times does `trim` run?
2. How many times does `align` run?
3. What does the `trimmed` channel contain?
4. Can `align` start before all `trim` jobs finish?

<details>
<summary>Click for answers</summary>

**Visual flow**:
```
reads channel: [file1.fastq] [file2.fastq] [file3.fastq]
                     ↓              ↓              ↓
                 [ trim ]       [ trim ]       [ trim ]    ← 3 parallel runs
                     ↓              ↓              ↓
trimmed channel: [trimmed.fq] [trimmed.fq] [trimmed.fq]
                     ↓              ↓              ↓
ref channel:    [genome.fa]    [genome.fa]    [genome.fa]  ← Same ref
                     ↓              ↓              ↓
                 [ align]       [ align]       [ align]    ← 3 parallel runs
                     ↓              ↓              ↓
aligned channel: [aligned.bam][aligned.bam][aligned.bam]
```

**Answers**:

1. **3 times** - Once for each input FASTQ file, potentially in parallel

2. **3 times** - Once for each trimmed file from the previous step

3. **The trimmed channel contains**: 3 trimmed.fq files (one from each trim run). These flow into the align process.

4. **Yes!** ✅ As soon as the first `trim` job finishes, the first `align` job can start. They don't all have to finish. This is **pipelining**:
   ```
   Time 0:  trim(file1)  starts
   Time 5:  trim(file1)  finishes → align(trim1) starts
   Time 10: trim(file2)  finishes → align(trim2) starts
   Time 15: align(trim1) finishes
   Time 20: trim(file3)  finishes → align(trim3) starts
   ```

This pipelining maximizes resource usage!
</details>

---

## 📝 Key Takeaways

Before moving to Day 5, ensure you understand:

✅ **Channels are streams**, not variables holding data  
✅ **Queue channels** are consumed once, emit multiple values  
✅ **Value channels** are reusable, emit one value  
✅ **Channels enable automatic parallelization** of processes  
✅ **Multiple items in a channel = multiple process runs**  
✅ **Reference data uses value channels** to be reused  
✅ **Sample data uses queue channels** to be processed  

### The Mental Model: Channel as Conveyor Belt

```
Channel → Process → Output Channel → Next Process → ...

[item1]─┐
[item2]─┼→ [Process] → [result1]─┐
[item3]─┘                [result2]─┼→ [Next Process]
                         [result3]─┘
```

Each item triggers a process run, results flow to the next stage.

### Common Channel Patterns

**Pattern 1: Files to process**
```groovy
files = Channel.fromPath("data/*.fastq")
```

**Pattern 2: Paired files**
```groovy
pairs = Channel.fromFilePairs("data/*_{1,2}.fastq")
```

**Pattern 3: Reference data**
```groovy
reference = Channel.value("genome.fa")
```

**Pattern 4: List of values**
```groovy
params = Channel.from(1, 2, 3, 4, 5)
```

**Pattern 5: From parameters**
```groovy
threads = Channel.value(params.cpus)
```

---

## 🎯 Ready for Day 5?

Tomorrow, you'll write your **first complete workflow** by connecting processes with channels! You'll see everything come together.

### Quick Preview

```groovy
// Tomorrow you'll write workflows like this:
workflow {
    // Create channels
    samples = Channel.fromPath("*.fastq")
    reference = Channel.value("genome.fa")
    
    // Chain processes
    qc_results = FASTQC(samples)
    aligned = ALIGN(samples, reference)
    variants = CALL_VARIANTS(aligned, reference)
    
    // Your first complete pipeline!
}
```

### Channel Quick Reference

```groovy
// Creating channels
Channel.fromPath("*.fastq")              // Files
Channel.fromFilePairs("*_{1,2}.fq")     // Paired files
Channel.from(1, 2, 3)                    // Values
Channel.value("reference.fa")            // Single reusable value
Channel.of(1, 2, 3)                      // Modern syntax for values

// Channel types
Queue channel:  consumed once, multiple values
Value channel:  reusable, single value

// In processes
process example {
    input:
    path fastq           // From queue channel
    path ref             // From value channel (or queue)
    
    output:
    path "*.bam"         // Creates new queue channel
}
```

---

## ✅ Day 4 Completion Checklist

Before marking Day 4 complete, ensure you can:

- [ ] Explain what a channel is (stream, not variable)
- [ ] Distinguish queue channels from value channels
- [ ] Create channels from files with `fromPath`
- [ ] Create channels from values with `value` or `from`
- [ ] Understand how channels enable parallelization
- [ ] Know when to use queue vs value channels
- [ ] Connect channels to process inputs
- [ ] Predict how many times a process will run

**Self-Test**: Given this scenario, create the appropriate channels:
- 50 BAM files to analyze
- 1 reference genome (used by all)
- Quality thresholds to test: 20, 30, 40

<details>
<summary>Check your answer</summary>

```groovy
// 50 BAM files - queue channel (different data)
bam_files = Channel.fromPath("data/*.bam")

// Reference genome - value channel (reusable)
reference = Channel.value("refs/genome.fa")

// Thresholds - queue channel (different values to test)
thresholds = Channel.from(20, 30, 40)

// Usage
workflow {
    // Each BAM × each threshold = 150 process runs
    results = analyzeBAM(bam_files, reference, thresholds)
}
```

If you got this right, you understand channels! 🎉
</details>

**Completed Day 4?** Update your `PROGRESS.md`! Channels are the key to Nextflow's power! 🚀

**Your progress**: 4/28 days (14.3%) complete

---

*Tomorrow: Day 5 - Connecting Processes into Workflows*

**Great work! See you tomorrow! 🚀**