# Day 1: What Nextflow Actually Is

**Learning Time**: 30 minutes  
**Prerequisites**: Basic Python knowledge, familiarity with bioinformatics concepts  
**Goal**: Understand what Nextflow is, why it exists, and when to use it

---

## 📖 Introduction (5 minutes)

Welcome to your Nextflow journey! Today, you'll learn what Nextflow is and, more importantly, *why* it exists. As a Python programmer, you might wonder: "Why learn another tool when I can script everything in Python?"

The answer lies in understanding the difference between **data processing** (what Python excels at) and **workflow orchestration** (what Nextflow was built for).

### What You'll Learn Today

- What Nextflow actually does and why bioinformaticians use it
- The fundamental problems Nextflow solves
- How Nextflow differs from writing Python scripts
- When to use Nextflow vs Python
- Real-world bioinformatics scenarios where Nextflow shines

---

## 🎯 What Is Nextflow? (8 minutes)

### The Simple Definition

**Nextflow is a workflow orchestration engine designed for data-intensive computational pipelines, particularly in bioinformatics.**

Let's break this down:

- **Workflow**: A series of computational steps (running tools like BWA, GATK, BLAST)
- **Orchestration**: Coordinating when and how these steps run, managing their inputs/outputs
- **Data-intensive**: Handling hundreds or thousands of samples efficiently
- **Engine**: The system that executes your pipeline definition

### The Problem Nextflow Solves

Imagine you're analyzing 100 genomic samples. Each sample needs:
1. Quality control (FastQC)
2. Alignment (BWA)
3. Variant calling (GATK)
4. Annotation (SnpEff)

**The Python Approach**:
```python
import subprocess
import glob

# Process each sample sequentially
for sample in glob.glob("samples/*.fastq"):
    # Step 1: QC
    subprocess.run(["fastqc", sample])
    
    # Step 2: Align
    subprocess.run(["bwa", "mem", "ref.fa", sample, "-o", f"{sample}.bam"])
    
    # Step 3: Variant calling
    subprocess.run(["gatk", "HaplotypeCaller", "-I", f"{sample}.bam"])
    
    # Step 4: Annotation
    subprocess.run(["snpeff", f"{sample}.vcf"])
```

**Problems with this approach**:
1. **Sequential execution**: Sample 2 waits for Sample 1 to finish completely
2. **No resumability**: If it crashes at sample 50, you restart from sample 1
3. **Manual parallelization**: You'd need to write threading/multiprocessing code
4. **Environment issues**: Different tools need different dependencies
5. **Not portable**: Might work on your laptop but not on the cluster
6. **Manual resource management**: You specify CPUs/memory manually for each run

**The Nextflow Approach**:
```groovy
// Define what each step does (we'll learn this syntax soon)
process fastqc {
    input: path(sample)
    output: path("*.html")
    
    "fastqc ${sample}"
}

process align {
    input: path(sample)
    output: path("*.bam")
    
    "bwa mem ref.fa ${sample} > ${sample}.bam"
}

process callVariants {
    input: path(bam)
    output: path("*.vcf")
    
    "gatk HaplotypeCaller -I ${bam}"
}

// Define how data flows through the workflow
workflow {
    samples = Channel.fromPath("samples/*.fastq")
    qc_results = fastqc(samples)
    aligned = align(samples)
    variants = callVariants(aligned)
}
```

**What Nextflow automatically handles**:
1. ✅ **Parallel execution**: All 100 samples run FastQC simultaneously (respecting system limits)
2. ✅ **Resumability**: `-resume` flag continues from where it failed
3. ✅ **Automatic parallelization**: No threading code needed
4. ✅ **Container support**: Each process can specify its own Docker/Singularity container
5. ✅ **Portability**: Same code runs on laptop, cluster, or cloud
6. ✅ **Resource management**: Declaratively specify CPU/memory requirements

### The Core Insight

**Python is great for *what* to do with data.**  
**Nextflow is great for *orchestrating* the doing.**

Think of it this way:
- **Python** = The chef who knows how to cook each dish
- **Nextflow** = The kitchen manager who coordinates multiple chefs, ensures dishes come out in order, and handles problems when a chef is sick

---

## 🔑 Key Concepts with Examples (10 minutes)

### 1. Declarative vs Imperative Programming

**Imperative (Python)**: You tell the computer *how* to do something step-by-step
```python
# You control the exact flow
results = []
for file in files:
    result = process(file)
    results.append(result)
```

**Declarative (Nextflow)**: You describe *what* you want to happen
```groovy
// Nextflow figures out how to execute efficiently
workflow {
    files = Channel.fromPath("*.txt")
    process_files(files)  // Automatically parallelized
}
```

### 2. Data Flow vs Control Flow

**Control Flow (Traditional Programming)**:
```python
# You explicitly control when things happen
step1_output = step1(input_data)
if step1_output.is_valid():
    step2_output = step2(step1_output)
    step3(step2_output)
```

**Data Flow (Nextflow)**:
```groovy
// Data flowing through processes triggers execution
workflow {
    input_ch = Channel.fromPath("data.txt")
    result1 = step1(input_ch)
    result2 = step2(result1)  // Automatically waits for step1
    step3(result2)            // Automatically waits for step2
}
```

The difference: In Nextflow, processes execute *when their input data is available*, not when you explicitly call them.

### 3. Processes as Isolated Units

**Python Function**:
```python
def align_reads(fastq_file):
    # Runs in the same Python environment
    # Shares memory space
    # Uses current working directory
    subprocess.run(["bwa", "mem", "ref.fa", fastq_file])
```

**Nextflow Process**:
```groovy
process alignReads {
    // Runs in isolated directory (work/xx/yyyy...)
    // Has its own compute resources
    // Can run in its own container
    // Can run on a different machine
    
    input:
    path fastq
    
    output:
    path "*.bam"
    
    script:
    """
    bwa mem ref.fa ${fastq} > aligned.bam
    """
}
```

Each Nextflow process runs in isolation, making it:
- **Reproducible**: Same inputs always produce same outputs
- **Portable**: Can run anywhere (local, cluster, cloud)
- **Parallelizable**: Can run many instances simultaneously
- **Resumable**: Failed processes can be retried independently

### 4. Channels: The Data Highways

**Python List (All in Memory)**:
```python
# Load everything at once
files = glob.glob("*.fastq")  # Could be 10,000 files!
for f in files:
    process(f)  # Sequential
```

**Nextflow Channel (Streaming)**:
```groovy
// Emits items as a stream
Channel
    .fromPath("*.fastq")
    .set { fastq_ch }

// Each file flows through and triggers parallel processing
process_fastq(fastq_ch)  // Automatically parallel for all files
```

Think of channels as conveyor belts in a factory: items (files, data) move along the belt, and workers (processes) handle them as they arrive.

### 5. Automatic Parallelization Example

Let's see a concrete comparison:

**Python - Manual Parallelization**:
```python
from concurrent.futures import ProcessPoolExecutor
import subprocess

def process_sample(sample):
    # Run quality control
    subprocess.run(["fastqc", sample])
    return sample

# Manually set up parallel processing
samples = glob.glob("*.fastq")
with ProcessPoolExecutor(max_workers=4) as executor:
    results = list(executor.map(process_sample, samples))

# Now you need to chain the next step...
# And handle errors, retries, monitoring...
```

**Nextflow - Automatic Parallelization**:
```groovy
process qualityControl {
    input:
    path sample
    
    output:
    path "*.html"
    
    script:
    """
    fastqc ${sample}
    """
}

workflow {
    samples = Channel.fromPath("*.fastq")
    qualityControl(samples)  // That's it! Automatic parallelization
}
```

Nextflow automatically:
- Distributes work across available CPUs
- Respects system resource limits
- Handles task scheduling
- Manages task dependencies
- Provides monitoring and logging

---

## 🔗 Python Connection: When to Use Each (3 minutes)

Understanding when to use Python vs Nextflow is crucial:

### Use Python When:
- ✅ Writing data processing algorithms (parsing, transforming, analyzing)
- ✅ Performing statistical analysis
- ✅ Creating visualizations
- ✅ Prototyping quick analysis scripts
- ✅ Working with single samples or small datasets
- ✅ Building custom analysis tools

**Example**: Computing GC content from a FASTA file
```python
from Bio import SeqIO

def calculate_gc_content(fasta_file):
    gc_count = 0
    total_count = 0
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq).upper()
        gc_count += sequence.count('G') + sequence.count('C')
        total_count += len(sequence)
    
    return (gc_count / total_count) * 100

# Perfect use of Python!
gc_percent = calculate_gc_content("genome.fasta")
print(f"GC content: {gc_percent:.2f}%")
```

### Use Nextflow When:
- ✅ Orchestrating multiple bioinformatics tools
- ✅ Processing many samples (10s to 1000s)
- ✅ Building reproducible pipelines
- ✅ Need to scale from laptop to cluster to cloud
- ✅ Want automatic parallelization and resumability
- ✅ Need to coordinate different software dependencies

**Example**: Processing 100 genomes through multiple tools
```groovy
// Perfect use of Nextflow!
workflow {
    genomes = Channel.fromPath("genomes/*.fasta")
    
    // Each step automatically parallelized across all genomes
    qc_results = qualityCheck(genomes)
    annotations = annotateGenes(qc_results)
    comparisons = compareGenomes(annotations)
    reports = generateReport(comparisons)
}
```

### The Best Approach: Use Both Together!

```groovy
// Nextflow orchestrates the pipeline
process customAnalysis {
    input:
    path data_file
    
    output:
    path "results.csv"
    
    script:
    """
    # Call your Python script for actual analysis
    python ${projectDir}/scripts/analyze_data.py ${data_file} > results.csv
    """
}
```

You can (and should) use Python scripts within Nextflow processes for complex data transformations, while Nextflow handles the orchestration.

---

## 💻 Hands-On Exercises (10 minutes)

### Exercise 1: Identify the Right Tool (3 minutes)

For each scenario, decide: **Python** or **Nextflow** or **Both**?

**Scenario A**: You need to parse a single VCF file and count the number of SNPs vs INDELs.

<details>
<summary>Click for answer</summary>

**Answer**: Python

This is a single-file data processing task. Python is perfect for this.

```python
def count_variants(vcf_file):
    snps = 0
    indels = 0
    
    with open(vcf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            ref = fields[3]
            alt = fields[4]
            
            if len(ref) == 1 and len(alt) == 1:
                snps += 1
            else:
                indels += 1
    
    return snps, indels
```
</details>

**Scenario B**: You have 500 RNA-seq samples that need to go through: quality control → trimming → alignment → quantification → differential expression.

<details>
<summary>Click for answer</summary>

**Answer**: Nextflow (orchestration) + Python (analysis scripts)

This is a multi-step pipeline with many samples. Nextflow orchestrates, Python can handle the differential expression analysis.

```groovy
workflow {
    samples = Channel.fromPath("samples/*.fastq")
    
    qc = FASTQC(samples)
    trimmed = TRIM(samples)
    aligned = ALIGN(trimmed)
    counts = QUANTIFY(aligned)
    
    // Call Python script for DE analysis
    DE_ANALYSIS(counts.collect())
}
```
</details>

**Scenario C**: You want to download 10 genome assemblies from NCBI and run a custom similarity analysis you wrote in Python.

<details>
<summary>Click for answer</summary>

**Answer**: Both

Nextflow orchestrates the downloads and coordinates running your Python script on each genome.

```groovy
process downloadGenome {
    output:
    path "*.fasta"
    
    script:
    """
    ncbi-genome-download ...
    """
}

process analyzeSimilarity {
    input:
    path genome
    
    output:
    path "results.txt"
    
    script:
    """
    python ${projectDir}/analyze.py ${genome} > results.txt
    """
}

workflow {
    genomes = downloadGenome()
    results = analyzeSimilarity(genomes)
}
```
</details>

### Exercise 2: Understanding Automatic Parallelization (4 minutes)

Look at these two code snippets and answer the questions:

**Python Version**:
```python
import time

samples = ["sample1.fastq", "sample2.fastq", "sample3.fastq"]

start = time.time()
for sample in samples:
    print(f"Processing {sample}...")
    time.sleep(2)  # Simulates 2 seconds of work
    print(f"Finished {sample}")
end = time.time()

print(f"Total time: {end - start:.1f} seconds")
```

**Questions**:
1. How long will this take to run?
2. If you had 100 samples, how long would it take?
3. How would you make it parallel?

<details>
<summary>Click for answers</summary>

**Answers**:
1. **6 seconds** (3 samples × 2 seconds each, sequential)
2. **200 seconds** (100 samples × 2 seconds each, sequential)
3. You'd need to use `concurrent.futures`, `multiprocessing`, or `threading`:

```python
from concurrent.futures import ThreadPoolExecutor
import time

def process_sample(sample):
    print(f"Processing {sample}...")
    time.sleep(2)
    print(f"Finished {sample}")
    return sample

samples = ["sample1.fastq", "sample2.fastq", "sample3.fastq"]

start = time.time()
with ThreadPoolExecutor(max_workers=3) as executor:
    executor.map(process_sample, samples)
end = time.time()

print(f"Total time: {end - start:.1f} seconds")  # ~2 seconds if 3 workers
```
</details>

**Nextflow Version**:
```groovy
process analyzeSample {
    input:
    val sample
    
    script:
    """
    echo "Processing ${sample}..."
    sleep 2
    echo "Finished ${sample}"
    """
}

workflow {
    samples = Channel.of("sample1.fastq", "sample2.fastq", "sample3.fastq")
    analyzeSample(samples)
}
```

**Questions**:
1. How long will this take to run (assuming you have 3+ CPU cores)?
2. If you had 100 samples (and 8 CPU cores), approximately how long?
3. What did you have to do to make it parallel?

<details>
<summary>Click for answers</summary>

**Answers**:
1. **~2 seconds** - All three run simultaneously (if you have 3+ cores)
2. **~25 seconds** - 100 samples ÷ 8 cores = 13 batches, each taking 2 seconds
3. **Nothing!** Nextflow parallelizes automatically based on:
   - Available CPU cores
   - System resources
   - Configuration settings

You just define WHAT to do, Nextflow handles HOW to parallelize it.
</details>

### Exercise 3: Recognizing Workflow Patterns (3 minutes)

Read this Python script and identify what makes it a good candidate for Nextflow:

```python
import subprocess
import os

samples = [
    "patient001_tumor",
    "patient001_normal", 
    "patient002_tumor",
    "patient002_normal",
    # ... 100 more samples
]

reference_genome = "hg38.fa"
known_variants = "dbsnp.vcf"

for sample in samples:
    # Step 1: Quality check (takes ~5 min per sample)
    print(f"Running FastQC on {sample}...")
    subprocess.run(["fastqc", f"{sample}.fastq"])
    
    # Step 2: Alignment (takes ~30 min per sample)
    print(f"Aligning {sample}...")
    subprocess.run([
        "bwa", "mem", reference_genome, 
        f"{sample}.fastq", "-o", f"{sample}.sam"
    ])
    
    # Step 3: Convert to BAM (takes ~5 min per sample)
    print(f"Converting {sample} to BAM...")
    subprocess.run([
        "samtools", "view", "-b", 
        f"{sample}.sam", "-o", f"{sample}.bam"
    ])
    
    # Step 4: Call variants (takes ~20 min per sample)
    print(f"Calling variants for {sample}...")
    subprocess.run([
        "gatk", "HaplotypeCaller",
        "-R", reference_genome,
        "-I", f"{sample}.bam",
        "-O", f"{sample}.vcf"
    ])
    
    print(f"Finished {sample}!")

print("All samples processed!")
```

**Questions**:
1. How long will this take for 100 samples?
2. What happens if it crashes at sample 50?
3. What problems might occur with tool dependencies?
4. How would you run this on a computing cluster?
5. List 3 specific benefits Nextflow would provide here.

<details>
<summary>Click for answers</summary>

**Answers**:

1. **~100 hours (4+ days)** 
   - 60 minutes per sample × 100 samples = 6,000 minutes sequential

2. **You start over from sample 1**
   - No checkpointing or resumability
   - All work on samples 1-49 is wasted time

3. **Dependency conflicts**:
   - FastQC might need Java 8
   - GATK might need Java 11
   - BWA and samtools might need different library versions
   - Everything must coexist in one environment

4. **Significant refactoring needed**:
   - Write job submission scripts (SLURM/PBS)
   - Manually split samples across jobs
   - Handle job dependencies
   - Collect results from distributed locations
   - Monitor job status manually

5. **Three Nextflow benefits**:
   
   **Benefit 1: Automatic Parallelization**
   - All 100 samples would process simultaneously (limited by cluster resources)
   - ~1 hour total instead of 100 hours

   **Benefit 2: Resumability**
   - If it crashes at sample 50, `-resume` continues from there
   - Completed work is never redone

   **Benefit 3: Container Isolation**
   - Each tool runs in its own container with exact dependencies
   - No version conflicts
   - Reproducible across any system

The Nextflow version would look like:
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    
    qc = FASTQC(samples)
    aligned = BWA_ALIGN(samples, ref_genome)
    bam = SAM_TO_BAM(aligned)
    variants = CALL_VARIANTS(bam, ref_genome)
}
```

And run with:
```bash
nextflow run pipeline.nf -resume  # Automatically parallel & resumable!
```
</details>

---

## 🤔 Reflection Activity (4 minutes)

Take a moment to think about your own work and answer these questions. Write your thoughts in your progress log.

### Question 1: Your Current Workflow
Think about a bioinformatics analysis you've done (or would like to do) that involved multiple steps or multiple samples.

**Describe it briefly**:
- What were the steps?
- How many samples?
- How did you coordinate the steps?

**Reflection prompt**: 
- What was frustrating about managing this workflow?
- What would you change if you could?

### Question 2: Python vs Nextflow
Look at your answer above.

**Ask yourself**:
- Which parts were true "data processing" (good for Python)?
- Which parts were "orchestration" (good for Nextflow)?
- Where did you spend most of your time? (Coding analysis vs managing workflow)

### Question 3: Learning Goals
Based on what you learned today:

**What are you most excited to learn in Nextflow?**
- [ ] Automatic parallelization
- [ ] Resumability (never restart from scratch)
- [ ] Portability (same code everywhere)
- [ ] Container integration
- [ ] Scaling to large datasets
- [ ] Something else: ________________

**What concerns do you have about learning Nextflow?**
- [ ] Learning Groovy syntax
- [ ] Understanding channels
- [ ] Getting it set up
- [ ] Knowing when to use it vs Python
- [ ] Other: ________________

---

## 📝 Key Takeaways

Before moving to Day 2, make sure you understand:

✅ **Nextflow is an orchestration tool**, not a replacement for Python  
✅ **Nextflow excels at coordinating multiple tools and samples** automatically  
✅ **Automatic parallelization** means you don't write threading code  
✅ **Resumability** means failed workflows continue where they stopped  
✅ **Portability** means the same code runs on laptop, cluster, or cloud  
✅ **Containers** solve the "works on my machine" problem  
✅ **Python and Nextflow work together**: Python processes, Nextflow orchestrates

### The Mental Model

Think of bioinformatics workflows like a restaurant kitchen:

- **Python** = The chefs who prepare each dish (data processing)
- **Nextflow** = The kitchen manager who coordinates orders, assigns chefs, handles problems (orchestration)
- **Processes** = Individual cooking stations (pasta, grill, salad)
- **Channels** = The order tickets flowing through the kitchen
- **Containers** = Each station having its own specialized equipment

You wouldn't ask a single chef to prepare 100 orders sequentially, and you wouldn't have the kitchen manager cook the food. Each has its role.

---

## 🎯 Ready for Day 2?

Tomorrow, you'll learn **Groovy essentials**—just enough to read and write Nextflow code comfortably. You'll see it's not as different from Python as you might think!

### Quick Prep for Tomorrow
Take 2 minutes to scan these Groovy vs Python examples:

**Python**:
```python
name = "Alice"
message = f"Hello, {name}!"
numbers = [x * 2 for x in range(5)]
```

**Groovy**:
```groovy
name = "Alice"
message = "Hello, ${name}!"
numbers = (0..4).collect { it * 2 }
```

See? Not so different! Tomorrow we'll make you comfortable with these patterns.

---

## ✅ Day 1 Completion Checklist

Before marking Day 1 complete, ensure you can:

- [ ] Explain what Nextflow is in one sentence
- [ ] Describe the difference between orchestration and processing
- [ ] List three problems Nextflow solves
- [ ] Explain when to use Python vs Nextflow
- [ ] Understand what automatic parallelization means
- [ ] Recognize a workflow that would benefit from Nextflow

**Completed Day 1?** Update your `PROGRESS.md` and celebrate! You've taken the first step toward mastering workflow orchestration. 🎉

**Time to Day 2**: 24 hours  
**Your progress**: 1/28 days (3.6%) complete

---

*Tomorrow: Day 2 - Groovy Essentials for Nextflow*

**See you tomorrow! 🚀**