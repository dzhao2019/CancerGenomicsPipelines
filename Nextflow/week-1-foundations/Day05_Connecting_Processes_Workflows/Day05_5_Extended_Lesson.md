# Day 5: Connecting Processes into Workflows

**Learning Time**: 30 minutes  
**Prerequisites**: Days 1-4 completed, understanding of processes and channels  
**Goal**: Write complete, multi-step workflows that orchestrate multiple processes

---

## 📖 Introduction (3 minutes)

Welcome to Day 5! This is where everything comes together. You've learned to write processes (Day 3) and work with channels (Day 4). Today, you'll **connect them into complete workflows** that solve real bioinformatics problems.

**What makes today exciting**: By the end of this session, you'll have written complete, working pipelines that chain multiple processes together—the kind of workflows you'll use in real research.

### What You'll Learn Today

- The workflow block syntax and structure
- How to chain processes together via channels
- Data flow patterns (linear, branching, merging)
- How process outputs automatically become channels
- Building complete multi-step pipelines
- The declarative workflow approach

### The Big Picture

**Individual processes** are like LEGO blocks:
```groovy
process FASTQC { ... }
process TRIM { ... }
process ALIGN { ... }
```

**Workflows** are how you assemble them:
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    qc = FASTQC(samples)
    trimmed = TRIM(samples)
    aligned = ALIGN(trimmed, reference)
}
```

Think of workflows as the **assembly instructions** that turn individual blocks into a complete structure.

---

## 🎯 Key Concepts with Examples (12 minutes)

### 1. The Workflow Block (2 minutes)

The `workflow` block is where you orchestrate your pipeline:

```groovy
workflow {
    // 1. Create input channels
    // 2. Call processes
    // 3. Connect outputs to inputs
    // 4. Define data flow
}
```

**Basic Structure**:
```groovy
// Define processes first
process stepOne {
    input: path(input_file)
    output: path("output.txt")
    script: "tool1 ${input_file} > output.txt"
}

process stepTwo {
    input: path(input_file)
    output: path("final.txt")
    script: "tool2 ${input_file} > final.txt"
}

// Then orchestrate in workflow block
workflow {
    // Create input channel
    input_ch = Channel.fromPath("data/*.txt")
    
    // Run first process
    result1 = stepOne(input_ch)
    
    // Run second process using first's output
    result2 = stepTwo(result1)
}
```

**Key Insights**:
- Processes are defined **outside** the workflow block
- The workflow block **calls** processes and **connects** their data
- Process outputs automatically become channels
- You can store results in variables for clarity

### 2. How Processes Connect: Output → Input (3 minutes)

This is the fundamental pattern you'll use constantly:

**Process outputs become channels automatically**:
```groovy
process createFile {
    output:
    path "data.txt"
    
    script:
    """
    echo "Hello" > data.txt
    """
}

workflow {
    // When you call a process, it returns a channel
    output_channel = createFile()
    
    // This channel can feed into another process
    nextProcess(output_channel)
}
```

**Chaining Processes**:
```groovy
process trim {
    input: path(fastq)
    output: path("trimmed.fq")
    script: "trim_tool ${fastq} > trimmed.fq"
}

process align {
    input: path(trimmed)
    output: path("aligned.bam")
    script: "bwa mem ref.fa ${trimmed} > aligned.bam"
}

process sort {
    input: path(bam)
    output: path("sorted.bam")
    script: "samtools sort ${bam} > sorted.bam"
}

workflow {
    reads = Channel.fromPath("*.fastq")
    
    // Chain processes together
    trimmed = trim(reads)
    aligned = align(trimmed)
    sorted = sort(aligned)
    
    // Data flows: reads → trim → align → sort
}
```

**What happens**:
1. `reads` channel emits FASTQ files
2. Each FASTQ goes through `trim`, producing `trimmed` channel
3. Each trimmed file goes through `align`, producing `aligned` channel
4. Each BAM goes through `sort`, producing `sorted` channel
5. All steps can pipeline (next step starts before previous finishes all items)

**Inline Style** (more compact):
```groovy
workflow {
    reads = Channel.fromPath("*.fastq")
    
    // Chain without intermediate variables
    sorted = sort(align(trim(reads)))
    
    // Same as above, reads from inside out:
    // reads → trim() → align() → sort()
}
```

**Pipe Operator Style** (most readable):
```groovy
workflow {
    Channel.fromPath("*.fastq")
        | trim
        | align
        | sort
}
```

All three styles do the same thing—choose what's most readable for your workflow!

### 3. Multiple Inputs in Workflows (2 minutes)

Many processes need multiple inputs—some from previous processes, some from reference data:

```groovy
process align {
    input:
    tuple val(sample_id), path(reads)
    path reference  // From a value channel
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    script:
    """
    bwa mem ${reference} ${reads} > ${sample_id}.bam
    """
}

process callVariants {
    input:
    tuple val(sample_id), path(bam)
    path reference
    path known_sites
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf")
    
    script:
    """
    gatk HaplotypeCaller \\
        -R ${reference} \\
        -I ${bam} \\
        --dbsnp ${known_sites} \\
        -O ${sample_id}.vcf
    """
}

workflow {
    // Queue channel: different per sample
    reads = Channel.fromFilePairs("data/*_{R1,R2}.fastq")
    
    // Value channels: same for all samples
    reference = Channel.value("refs/genome.fa")
    known_sites = Channel.value("refs/dbsnp.vcf")
    
    // Chain processes with multiple inputs
    bam_ch = align(reads, reference)
    vcf_ch = callVariants(bam_ch, reference, known_sites)
}
```

**Pattern**: 
- **Queue channels** for sample-specific data (flows through pipeline)
- **Value channels** for shared reference data (used by all samples)

### 4. Branching and Merging Workflows (2 minutes)

Sometimes data needs to split into multiple paths or merge from multiple sources:

**Branching** (one input, multiple processes):
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    
    // Same input goes to multiple processes
    qc_results = FASTQC(samples)
    kmer_analysis = KMER_COUNT(samples)
    read_stats = READ_STATS(samples)
    
    // All three processes run in parallel on the same input!
}
```

⚠️ **Important**: This works because each process call consumes from the channel **simultaneously**, not sequentially.

**Merging** (multiple inputs, one process):
```groovy
process combineResults {
    input:
    path "qc/*"
    path "stats/*"
    
    output:
    path "combined_report.html"
    
    script:
    """
    combine_reports.sh qc/ stats/ > combined_report.html
    """
}

workflow {
    samples = Channel.fromPath("*.fastq")
    
    qc = FASTQC(samples)
    stats = READ_STATS(samples)
    
    // Collect all outputs before combining
    combined = combineResults(
        qc.collect(),
        stats.collect()
    )
}
```

**The `.collect()` operator**: Waits for all items from a channel, then emits them as a single collection.

**Advanced Branching with Tuple Outputs**:
```groovy
process splitAnalysis {
    input:
    path data
    
    output:
    path "results.txt", emit: results
    path "report.html", emit: report
    path "summary.csv", emit: summary
    
    script:
    """
    analyze.sh ${data}
    """
}

workflow {
    data = Channel.fromPath("*.txt")
    splitAnalysis(data)
    
    // Access specific outputs
    processResults(splitAnalysis.out.results)
    publishReport(splitAnalysis.out.report)
    createSummary(splitAnalysis.out.summary)
}
```

### 5. A Complete Real-World Workflow (3 minutes)

Let's build a complete RNA-seq quality control workflow:

```groovy
// Process 1: Quality Check
process FASTQC {
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_fastqc.html")
    
    script:
    """
    fastqc ${reads} -o .
    mv *_fastqc.html ${sample_id}_fastqc.html
    """
}

// Process 2: Trim Low Quality Bases
process TRIM {
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq")
    
    script:
    """
    trimmomatic SE ${reads} ${sample_id}_trimmed.fastq \\
        TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Process 3: Post-Trim QC
process FASTQC_TRIMMED {
    input:
    tuple val(sample_id), path(trimmed)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed_fastqc.html")
    
    script:
    """
    fastqc ${trimmed} -o .
    mv *_fastqc.html ${sample_id}_trimmed_fastqc.html
    """
}

// Process 4: Aggregate All QC Reports
process MULTIQC {
    publishDir "results/qc", mode: 'copy'
    
    input:
    path "*"
    
    output:
    path "multiqc_report.html"
    
    script:
    """
    multiqc .
    """
}

// Workflow: Orchestrate Everything
workflow {
    // Input: Paired sample ID with FASTQ files
    reads_ch = Channel.fromFilePairs("data/*_{R1,R2}.fastq")
        .map { id, files -> [id, files[0]] }  // Use R1 only for simplicity
    
    // Step 1: Initial QC
    raw_qc = FASTQC(reads_ch)
    
    // Step 2: Trim reads
    trimmed = TRIM(reads_ch)
    
    // Step 3: Post-trim QC
    trimmed_qc = FASTQC_TRIMMED(trimmed)
    
    // Step 4: Aggregate all QC reports
    all_qc = raw_qc
        .map { id, html -> html }  // Extract just the HTML files
        .mix(trimmed_qc.map { id, html -> html })
        .collect()
    
    MULTIQC(all_qc)
}
```

**What this workflow does**:
1. Takes paired-end FASTQ files
2. Runs FastQC on raw reads
3. Trims low-quality bases
4. Runs FastQC on trimmed reads
5. Aggregates all QC reports into one MultiQC report

**Data flow**:
```
reads_ch ──┬─→ FASTQC ─────────→ raw_qc ────┐
           │                                  │
           └─→ TRIM ──→ FASTQC_TRIMMED ──→ trimmed_qc ─→ mix → collect → MULTIQC
```

**New operators used**:
- `.map {}`: Transform channel items
- `.mix()`: Combine multiple channels
- `.collect()`: Gather all items into one

---

## 🔗 Python Connection: Scripts vs Workflows (3 minutes)

Let's compare how you'd orchestrate a multi-step pipeline:

### Python Script Approach

```python
import subprocess
import glob
import os

# Sequential pipeline script
def run_pipeline():
    samples = glob.glob("data/*.fastq")
    
    for sample in samples:
        sample_id = os.path.basename(sample).replace(".fastq", "")
        
        # Step 1: FastQC
        print(f"Running FastQC on {sample_id}...")
        subprocess.run(["fastqc", sample, "-o", "qc/"])
        
        # Step 2: Trimming
        print(f"Trimming {sample_id}...")
        trimmed = f"trimmed/{sample_id}_trimmed.fastq"
        subprocess.run([
            "trimmomatic", "SE", sample, trimmed,
            "TRAILING:3", "MINLEN:36"
        ])
        
        # Step 3: Post-trim QC
        print(f"Running FastQC on trimmed {sample_id}...")
        subprocess.run(["fastqc", trimmed, "-o", "qc/"])
    
    # Step 4: Aggregate
    print("Creating MultiQC report...")
    subprocess.run(["multiqc", "qc/", "-o", "results/"])

if __name__ == "__main__":
    run_pipeline()
```

**Issues**:
- ❌ Sequential processing (one sample at a time)
- ❌ Steps must complete in order (no pipelining)
- ❌ Hard to resume if it fails
- ❌ Manual error handling needed
- ❌ Hard to parallelize
- ❌ Difficult to track what succeeded/failed

### Nextflow Workflow Approach

```groovy
workflow {
    // Parallel processing automatically
    reads_ch = Channel.fromPath("data/*.fastq")
        .map { file -> [file.baseName, file] }
    
    // Pipeline data flow
    raw_qc = FASTQC(reads_ch)
    trimmed = TRIM(reads_ch)
    trimmed_qc = FASTQC_TRIMMED(trimmed)
    
    all_qc = raw_qc.mix(trimmed_qc).collect()
    MULTIQC(all_qc)
}
```

**Benefits**:
- ✅ Automatic parallelization
- ✅ Pipelined execution (steps overlap)
- ✅ Built-in resume with `-resume`
- ✅ Automatic error handling and retry
- ✅ Clear data flow
- ✅ Automatic tracking of all tasks

### Time Comparison

For 10 samples, 5 minutes per step:

| Approach | Total Time | Notes |
|----------|-----------|-------|
| Python sequential | 150 min (2.5 hrs) | 10 samples × 3 steps × 5 min |
| Nextflow (4 cores) | ~40 min | Parallel + pipelined |
| Nextflow (cluster) | ~15 min | Maximum parallelization |

### Code Comparison

| Aspect | Python | Nextflow |
|--------|--------|----------|
| **Lines of code** | ~30 lines | ~80 lines (more verbose) |
| **Parallelization** | Manual (complex) | Automatic |
| **Error recovery** | Manual checkpoints | `nextflow run -resume` |
| **Clarity** | Imperative steps | Declarative flow |
| **Scalability** | Single machine | Laptop to cloud |
| **Reproducibility** | Environment-dependent | Container-based |

**When to use each**:
- **Python**: Quick scripts, data analysis, single samples
- **Nextflow**: Multi-step pipelines, many samples, production workflows

---

## 💻 Hands-On Exercises (10 minutes)

### Exercise 1: Chain Two Processes (3 minutes)

**Task**: Create a simple two-step workflow.

**Given**:
```groovy
process countLines {
    input:
    path input_file
    
    output:
    path "line_count.txt"
    
    script:
    """
    wc -l ${input_file} > line_count.txt
    """
}

process summarize {
    input:
    path count_file
    
    output:
    path "summary.txt"
    
    script:
    """
    echo "Total files processed: \$(cat ${count_file} | wc -l)" > summary.txt
    """
}
```

**Your task**: Write a workflow that:
1. Creates a channel from files in `data/` directory
2. Counts lines in each file
3. Summarizes the results

<details>
<summary>Click for solution</summary>

```groovy
workflow {
    // Step 1: Create input channel
    files = Channel.fromPath("data/*.txt")
    
    // Step 2: Count lines in each file
    counts = countLines(files)
    
    // Step 3: Collect all counts and summarize
    summary = summarize(counts.collect())
    
    // Or as a chain:
    // Channel.fromPath("data/*.txt")
    //     | countLines
    //     | collect
    //     | summarize
}
```

**Explanation**:
- `files` channel emits each .txt file
- `countLines` runs once per file (parallel)
- `counts` channel contains all line_count.txt files
- `.collect()` waits for all counts, then emits as one collection
- `summarize` runs once with all collected counts

**Data flow**:
```
file1.txt ──→ countLines ──→ count1.txt ──┐
file2.txt ──→ countLines ──→ count2.txt ──┼─→ collect ──→ summarize
file3.txt ──→ countLines ──→ count3.txt ──┘
```
</details>

### Exercise 2: Workflow with Reference Data (3 minutes)

**Task**: Build an alignment workflow with reference genome.

**Given these processes**:
```groovy
process INDEX {
    input:
    path genome
    
    output:
    path "${genome}.*"
    
    script:
    """
    bwa index ${genome}
    """
}

process ALIGN {
    input:
    tuple val(sample_id), path(reads)
    path genome
    path "*"  // Index files
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    script:
    """
    bwa mem ${genome} ${reads} | samtools view -b > ${sample_id}.bam
    """
}
```

**Your task**: Write a workflow that:
1. Indexes the reference genome once
2. Aligns all samples using the indexed genome

**Hint**: The genome is used by all samples (value channel!)

<details>
<summary>Click for solution</summary>

```groovy
workflow {
    // Input samples (queue channel)
    samples = Channel.fromPath("data/*.fastq")
        .map { file -> [file.baseName, file] }
    
    // Reference genome (value channel)
    genome = Channel.value("refs/genome.fa")
    
    // Step 1: Index genome once
    index_files = INDEX(genome)
    
    // Step 2: Align all samples
    // Each sample gets the same genome and index files
    aligned = ALIGN(samples, genome, index_files)
}
```

**Explanation**:
- `genome` is a value channel → can be reused
- `INDEX` runs **once** (only one genome)
- `index_files` output is also reusable
- `ALIGN` runs once per sample, all use the same genome and index
- If you have 100 samples, INDEX runs 1×, ALIGN runs 100× (in parallel)

**Data flow**:
```
genome.fa ──→ INDEX ──→ index_files ──┐
                                      │
sample1 ──────────────────────────────┼──→ ALIGN ──→ sample1.bam
sample2 ──────────────────────────────┼──→ ALIGN ──→ sample2.bam
sample3 ──────────────────────────────┼──→ ALIGN ──→ sample3.bam
                                      │
                    (all use same index)
```
</details>

### Exercise 3: Complete QC Pipeline (4 minutes)

**Task**: Build a complete quality control workflow.

**Requirements**:
1. Input: FASTQ files from `samples/` directory
2. Run FastQC on each sample
3. Run a custom statistics script on each sample
4. Combine all QC reports with MultiQC
5. Publish final report to `results/` directory

**Given processes** (simplified):
```groovy
process FASTQC {
    input:
    tuple val(id), path(fastq)
    output:
    path "*.html"
    script:
    """
    fastqc ${fastq}
    """
}

process STATS {
    input:
    tuple val(id), path(fastq)
    output:
    path "${id}_stats.txt"
    script:
    """
    seqkit stats ${fastq} > ${id}_stats.txt
    """
}

process MULTIQC {
    publishDir "results", mode: 'copy'
    
    input:
    path "*"
    output:
    path "multiqc_report.html"
    script:
    """
    multiqc .
    """
}
```

**Write the complete workflow**:

<details>
<summary>Click for solution</summary>

```groovy
workflow {
    // Step 1: Create input channel with sample IDs
    samples = Channel.fromPath("samples/*.fastq")
        .map { file -> [file.baseName, file] }
    
    // Step 2: Run FastQC on all samples (parallel)
    fastqc_results = FASTQC(samples)
    
    // Step 3: Run stats on all samples (parallel)
    stats_results = STATS(samples)
    
    // Step 4: Combine all outputs for MultiQC
    all_qc = fastqc_results
        .mix(stats_results)
        .collect()
    
    // Step 5: Create aggregate report
    MULTIQC(all_qc)
}
```

**Alternative with pipe operator**:
```groovy
workflow {
    samples = Channel.fromPath("samples/*.fastq")
        .map { file -> [file.baseName, file] }
    
    // Branch into two parallel analyses
    fastqc_out = samples | FASTQC
    stats_out = samples | STATS
    
    // Merge and aggregate
    fastqc_out
        .mix(stats_out)
        .collect()
        | MULTIQC
}
```

**Explanation**:
- `samples` feeds both FASTQC and STATS simultaneously
- Both processes run in parallel across all samples
- `.mix()` combines the two output channels
- `.collect()` waits for all outputs, then emits as collection
- MULTIQC receives all reports at once

**Data flow visualization**:
```
sample1 ──┬──→ FASTQC ──→ sample1.html ──┐
sample2 ──┤                               ├──→ mix ──→ collect ──→ MULTIQC
sample3 ──┤                               │
          └──→ STATS ───→ sample1_stats ──┘
                          sample2_stats
                          sample3_stats
```

**Why this works**:
- Each sample independently flows through both processes
- Everything runs in parallel
- Results are collected at the end
- Clean, declarative workflow
</details>

---

## 🤔 Reflection Activity (4 minutes)

### Question 1: Understanding Data Flow

Consider this workflow:
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    
    step1_out = PROCESS_A(samples)
    step2_out = PROCESS_B(step1_out)
    step3_out = PROCESS_C(step1_out)  // Note: uses step1_out again!
}
```

**Questions**:
1. Is this valid? Will it work?
2. What happens to the `step1_out` channel?
3. How many times does each process run (if there are 5 FASTQ files)?

<details>
<summary>Click for answers</summary>

**Answers**:

1. **Is this valid?**
   
   ❌ **NO!** This will fail with a channel consumption error.

2. **What happens?**
   
   `step1_out` is a queue channel that gets consumed by `PROCESS_B`. When `PROCESS_C` tries to use it, the channel is already empty.
   
   **Error**: "Channel has been used as input by more than one process or operator"

3. **How to fix it?**

   **Option 1**: Use the output multiple times by splitting
   ```groovy
   workflow {
       samples = Channel.fromPath("*.fastq")
       
       step1_out = PROCESS_A(samples)
       
       // Duplicate the channel
       step1_out.multiMap { item ->
           branch_b: item
           branch_c: item
       }.set { branches }
       
       step2_out = PROCESS_B(branches.branch_b)
       step3_out = PROCESS_C(branches.branch_c)
   }
   ```
   
   **Option 2**: Call PROCESS_A twice (wasteful!)
   ```groovy
   workflow {
       samples = Channel.fromPath("*.fastq")
       
       step2_out = PROCESS_B(PROCESS_A(samples))
       step3_out = PROCESS_C(PROCESS_A(samples))  // Runs A again!
   }
   ```
   
   **Option 3**: Make processes sequential if possible
   ```groovy
   workflow {
       samples = Channel.fromPath("*.fastq")
       
       step1_out = PROCESS_A(samples)
       step2_out = PROCESS_B(step1_out)
       step3_out = PROCESS_C(step2_out)  // Uses B's output
   }
   ```

**Each process runs**: 5 times (once per FASTQ file), assuming successful splitting
</details>

### Question 2: Workflow Design

You need to build this pipeline:
1. Quality check raw reads
2. Trim adapters
3. Quality check trimmed reads
4. Align to reference
5. Call variants

**Design questions**:
- Which processes can run in parallel?
- Which must be sequential?
- Which need reference data (value channels)?
- Draw the data flow

<details>
<summary>Click for solution</summary>

**Parallel vs Sequential**:

**Can run in parallel**:
- QC on raw reads (step 1) - independent per sample
- All samples through the entire pipeline simultaneously

**Must be sequential** (per sample):
- Raw QC → Trim → Trimmed QC → Align → Variants (one sample's journey)

**Data flow**:
```
ref.fa ────────────────────────────┐
                                   │
sample1.fq ──→ QC1 ──→ TRIM ──→ QC2 ──→ ALIGN ──→ VARIANTS ──→ sample1.vcf
sample2.fq ──→ QC1 ──→ TRIM ──→ QC2 ──→ ALIGN ──→ VARIANTS ──→ sample2.vcf
sample3.fq ──→ QC1 ──→ TRIM ──→ QC2 ──→ ALIGN ──→ VARIANTS ──→ sample3.vcf
                                   │
                          (all use same ref)
```

**Workflow code**:
```groovy
workflow {
    // Queue: different per sample
    samples = Channel.fromPath("data/*.fastq")
        .map { file -> [file.baseName, file] }
    
    // Value: same for all samples
    reference = Channel.value("refs/genome.fa")
    known_sites = Channel.value("refs/dbsnp.vcf")
    
    // Linear pipeline per sample
    raw_qc = QC_RAW(samples)
    trimmed = TRIM(samples)
    trimmed_qc = QC_TRIMMED(trimmed)
    aligned = ALIGN(trimmed, reference)
    variants = CALL_VARIANTS(aligned, reference, known_sites)
}
```

**Parallelization**:
- All samples flow through simultaneously
- Within each sample: sequential steps
- Across samples: fully parallel
- If 100 samples and 100 cores: ~same time as 1 sample!
</details>

### Question 3: Collect or Not?

For each scenario, decide if you need `.collect()`:

**Scenario A**: Creating a single summary report from 100 sample analyses
**Scenario B**: Running variant calling on 100 BAM files
**Scenario C**: Generating a MultiQC report from all FastQC outputs
**Scenario D**: Aligning 100 FASTQ files to a reference

<details>
<summary>Click for answers</summary>

**Scenario A: YES** ✅ Need `.collect()`
```groovy
sample_results = ANALYZE(samples)
summary = CREATE_SUMMARY(sample_results.collect())
```
**Why**: The summary needs ALL results at once

**Scenario B: NO** ❌ Don't use `.collect()`
```groovy
bam_files = Channel.fromPath("*.bam")
variants = CALL_VARIANTS(bam_files, reference)
```
**Why**: Each BAM is processed independently

**Scenario C: YES** ✅ Need `.collect()`
```groovy
fastqc_results = FASTQC(samples)
multiqc_report = MULTIQC(fastqc_results.collect())
```
**Why**: MultiQC needs all QC files together

**Scenario D: NO** ❌ Don't use `.collect()`
```groovy
fastq_files = Channel.fromPath("*.fastq")
aligned = ALIGN(fastq_files, reference)
```
**Why**: Each FASTQ is aligned independently

**The Pattern**:
- **Aggregation/Summary** → Use `.collect()`
- **Independent processing** → No `.collect()`
</details>

---

## 📝 Key Takeaways

Before moving to Day 6, ensure you understand:

✅ **The workflow block orchestrates processes**  
✅ **Process outputs automatically become channels**  
✅ **Chaining: output of one process → input of next**  
✅ **Queue channels for sample data, value channels for references**  
✅ **`.collect()` for aggregation, not for independent processing**  
✅ **Data flow is declarative**—you describe what, not how  
✅ **Parallelization happens automatically** across samples  

### The Mental Model: Assembly Line

```
Raw        Quality      Trimming     Alignment    Variant
Materials  Check        Station      Station      Calling
  ↓          ↓             ↓            ↓            ↓
[sample] → [QC] ────→ [TRIM] ────→ [ALIGN] ───→ [VARIANTS]
[sample] → [QC] ────→ [TRIM] ────→ [ALIGN] ───→ [VARIANTS]
[sample] → [QC] ────→ [TRIM] ────→ [ALIGN] ───→ [VARIANTS]
           └─────────────────────────────────────────┘
                   All running simultaneously
```

Each sample moves through the assembly line, multiple samples at different stages concurrently.

### Common Workflow Patterns

**Pattern 1: Linear pipeline**
```groovy
workflow {
    input | step1 | step2 | step3
}
```

**Pattern 2: With reference data**
```groovy
workflow {
    samples = Channel.fromPath("*.fq")
    ref = Channel.value("genome.fa")
    samples | process(ref) | nextStep
}
```

**Pattern 3: Branching**
```groovy
workflow {
    samples = Channel.fromPath("*.fq")
    qc = FASTQC(samples)
    stats = STATS(samples)
    mix(qc, stats) | REPORT
}
```

**Pattern 4: Aggregation**
```groovy
workflow {
    results = ANALYZE(samples)
    summary = SUMMARIZE(results.collect())
}
```

---

## 🎯 Ready for Day 6?

Tomorrow, you'll learn how to **run workflows** using the Nextflow command line and understand what happens during execution!

### Quick Preview

```bash
# You'll learn to run your workflows:
nextflow run workflow.nf

# With parameters:
nextflow run workflow.nf --input data/ --output results/

# And resume failed runs:
nextflow run workflow.nf -resume

# With execution reports:
nextflow run workflow.nf -with-report report.html
```

### Workflow Quick Reference

```groovy
// Basic workflow structure
workflow {
    // 1. Create channels
    samples = Channel.fromPath("*.fastq")
    reference = Channel.value("genome.fa")
    
    // 2. Call processes (outputs become channels)
    qc = FASTQC(samples)
    aligned = ALIGN(samples, reference)
    
    // 3. Chain processes
    sorted = SORT(aligned)
    
    // 4. Aggregate if needed
    report = REPORT(qc.collect())
}

// Pipe operator style
workflow {
    Channel.fromPath("*.fastq")
        | QC
        | TRIM
        | ALIGN
}

// With branching
workflow {
    samples = Channel.fromPath("*.fastq")
    qc = QC(samples)
    stats = STATS(samples)
    qc.mix(stats).collect() | REPORT
}
```

---

## ✅ Day 5 Completion Checklist

Before marking Day 5 complete, ensure you can:

- [ ] Write a workflow block that orchestrates processes
- [ ] Chain processes together via their outputs
- [ ] Use value channels for reference data
- [ ] Use queue channels for sample data
- [ ] Understand when to use `.collect()`
- [ ] Read and understand workflow data flow
- [ ] Write a complete multi-step pipeline

**Self-Test**: Write a workflow for this scenario:
- FASTQ files in `data/`
- Reference genome at `refs/genome.fa`
- Pipeline: QC → Align → Sort → Index

<details>
<summary>Check your answer</summary>

```groovy
workflow {
    // Inputs
    samples = Channel.fromPath("data/*.fastq")
        .map { file -> [file.baseName, file] }
    reference = Channel.value("refs/genome.fa")
    
    // Pipeline
    qc = QC(samples)
    aligned = ALIGN(samples, reference)
    sorted = SORT(aligned)
    indexed = INDEX(sorted)
}

// Or with pipes:
workflow {
    samples = Channel.fromPath("data/*.fastq")
        .map { file -> [file.baseName, file] }
    reference = Channel.value("refs/genome.fa")
    
    samples
        | QC
        | { ALIGN(it, reference) }
        | SORT
        | INDEX
}
```

If you got this right, you can build workflows! 🎉
</details>

**Completed Day 5?** Update your `PROGRESS.md`! You can now build complete pipelines! 🚀

**Your progress**: 5/28 days (17.9%) complete

**🎊 Milestone achieved!** You've completed Week 1 foundations! You can now:
- Understand what Nextflow is and why it's valuable
- Read and write Groovy syntax
- Create Nextflow processes
- Work with channels
- Build complete workflows

Take a moment to celebrate—you've built a solid foundation! 🎉

---

*Tomorrow: Day 6 - Running Workflows and Understanding Execution*

**Excellent progress! See you tomorrow! 🚀**