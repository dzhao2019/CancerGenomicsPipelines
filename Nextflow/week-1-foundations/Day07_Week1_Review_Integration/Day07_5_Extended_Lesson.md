# Day 7: Week 1 Review and Integration Project

**Learning Time**: 30 minutes  
**Prerequisites**: Days 1-6 completed  
**Goal**: Consolidate Week 1 learning through a complete integration project

---

## 📖 Introduction (3 minutes)

Welcome to Day 7 - the final day of Week 1! 🎉

You've learned an incredible amount this week. Today is about **consolidation** - bringing everything together into one cohesive understanding. You'll review the key concepts from each day and then build a complete, production-style workflow that uses everything you've learned.

### Your Week 1 Journey

**Day 1**: What Nextflow Actually Is
- ✅ Understood orchestration vs processing
- ✅ Learned when to use Nextflow vs Python
- ✅ Grasped automatic parallelization concept

**Day 2**: Groovy Essentials
- ✅ Mastered string interpolation `${}`
- ✅ Learned closures and `it`
- ✅ Worked with Lists and Maps

**Day 3**: Your First Process
- ✅ Wrote processes with input/output/script
- ✅ Understood process isolation
- ✅ Used val, path, and tuple

**Day 4**: Understanding Channels
- ✅ Grasped channels as streams
- ✅ Learned queue vs value channels
- ✅ Created channels from various sources

**Day 5**: Connecting Processes into Workflows
- ✅ Wrote complete workflows
- ✅ Chained processes together
- ✅ Used operators like collect() and mix()

**Day 6**: Running and Understanding Execution
- ✅ Understood work directory structure
- ✅ Learned resume mechanism
- ✅ Debugged workflows effectively

### Today's Goals

1. **Quick review** of each day's core concept
2. **Build a complete integration project** using all concepts
3. **Self-assessment** to identify strengths and gaps
4. **Celebrate** your Week 1 achievement!
5. **Preview** Week 2 topics

---

## 🔄 Week 1 Concept Review (8 minutes)

Let's quickly review the essential concepts from each day. Test yourself - can you answer these without looking back?

### Day 1: The Core Understanding

**Key Question**: What's the fundamental difference between Nextflow and Python for bioinformatics?

<details>
<summary>Click to review</summary>

**Answer**: 
- **Python**: Data processing - great for algorithms, analysis, transformations
- **Nextflow**: Workflow orchestration - coordinates multiple tools, handles parallelization, manages execution

**Mental Model**: Python is the chef, Nextflow is the kitchen manager.

**When to use Nextflow**:
- Multiple bioinformatics tools to coordinate
- Many samples to process (10s to 1000s)
- Need automatic parallelization
- Want reproducibility and portability
- Building production pipelines
</details>

### Day 2: Essential Groovy Patterns

**Challenge**: Complete these Groovy snippets

```groovy
// 1. String interpolation with sample_id
filename = // sample_id is "patient001", create "patient001_results.bam"

// 2. Transform list to uppercase
samples = ["sample1", "sample2", "sample3"]
upper = samples. // → ["SAMPLE1", "SAMPLE2", "SAMPLE3"]

// 3. Filter numbers greater than 10
numbers = [5, 12, 8, 15, 3, 20]
big = numbers. // → [12, 15, 20]
```

<details>
<summary>Click for answers</summary>

```groovy
// 1. String interpolation
filename = "${sample_id}_results.bam"

// 2. Transform to uppercase
upper = samples.collect { it.toUpperCase() }

// 3. Filter
big = numbers.findAll { it > 10 }
```

**Key patterns to remember**:
- `"${variable}"` for interpolation
- `.collect {}` for transformation
- `.findAll {}` for filtering
- `it` for implicit parameter
</details>

### Day 3: Process Structure

**Challenge**: Write the skeleton of a process that aligns reads

<details>
<summary>Click for answer</summary>

```groovy
process alignReads {
    input:
    tuple val(sample_id), path(reads)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    script:
    """
    bwa mem ${reference} ${reads} | samtools view -b > ${sample_id}.bam
    """
}
```

**Key elements**:
- `input:` declares what data comes in
- `output:` declares what gets created
- `script:` contains shell commands with `${}` interpolation
- `tuple` keeps related data together
</details>

### Day 4: Channel Types

**Challenge**: For each scenario, choose queue or value channel

1. Reference genome used by all samples: ______
2. 100 FASTQ files to process: ______
3. List of quality thresholds to test: ______
4. Output directory parameter: ______

<details>
<summary>Click for answers</summary>

1. **Value channel** - reusable reference data
   ```groovy
   ref = Channel.value("genome.fa")
   ```

2. **Queue channel** - different files to process
   ```groovy
   fastqs = Channel.fromPath("*.fastq")
   ```

3. **Queue channel** - different values to test
   ```groovy
   thresholds = Channel.from(20, 30, 40)
   ```

4. **Value channel** - single config value
   ```groovy
   outdir = Channel.value(params.output)
   ```

**Remember**:
- Queue = consumed once, multiple values
- Value = reusable, single value
</details>

### Day 5: Workflow Patterns

**Challenge**: Complete this workflow structure

```groovy
workflow {
    // 1. Create channel from FASTQ files
    samples = 
    
    // 2. Create reference channel
    reference = 
    
    // 3. Run QC
    qc = 
    
    // 4. Align with reference
    aligned = 
    
    // 5. Collect all results
    all_results = 
}
```

<details>
<summary>Click for answer</summary>

```groovy
workflow {
    // 1. Create channel from FASTQ files
    samples = Channel.fromPath("data/*.fastq")
        .map { file -> [file.baseName, file] }
    
    // 2. Create reference channel
    reference = Channel.value("refs/genome.fa")
    
    // 3. Run QC
    qc = FASTQC(samples)
    
    // 4. Align with reference
    aligned = ALIGN(samples, reference)
    
    // 5. Collect all results
    all_results = aligned.collect()
}
```

**Key patterns**:
- Queue channels for sample data
- Value channels for shared references
- Process outputs become channels
- `.collect()` for aggregation
</details>

### Day 6: Execution Understanding

**Challenge**: True or False?

1. Each process task runs in the same directory: ____
2. The `-resume` flag reruns everything from scratch: ____
3. `.command.err` contains the error output: ____
4. Failed tasks are never cached: ____
5. Changing a parameter invalidates all cached tasks: ____

<details>
<summary>Click for answers</summary>

1. **FALSE** - Each task runs in its own `work/xx/yyyy/` directory
2. **FALSE** - Resume skips cached successful tasks
3. **TRUE** - Standard error goes to `.command.err`
4. **TRUE** - Failed tasks always retry on resume
5. **TRUE** - Parameters affect task hash

**Debug workflow**:
1. Find error in `.nextflow.log`
2. Navigate to `work/xx/yyyy/`
3. Check `.command.err`
4. Re-run with `bash .command.sh`
5. Fix and `nextflow run -resume`
</details>

---

## 🏗️ Integration Project: Complete Quality Control Pipeline (15 minutes)

Now let's build a **complete, production-style workflow** that integrates everything from Week 1.

### Project Overview

**Goal**: Create a comprehensive quality control pipeline for RNA-seq data

**Pipeline Steps**:
1. Initial quality check (FastQC)
2. Adapter trimming (Trimmomatic)
3. Post-trim quality check (FastQC)
4. Read statistics (custom script)
5. Aggregate report (MultiQC)

**What This Uses**:
- ✅ Multiple processes (Day 3)
- ✅ Queue and value channels (Day 4)
- ✅ Workflow orchestration (Day 5)
- ✅ Groovy syntax throughout (Day 2)
- ✅ Understanding of execution (Day 6)

### Complete Workflow Code

```groovy
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Week 1 Integration Project: RNA-seq Quality Control Pipeline
 * 
 * This workflow demonstrates all concepts from Week 1:
 * - Processes with proper input/output declarations
 * - Queue and value channels
 * - Workflow orchestration
 * - Data flow through multiple steps
 * - Aggregation and reporting
 */

// ============================================================================
// PARAMETERS
// ============================================================================

params.reads = "data/*_R{1,2}.fastq.gz"
params.outdir = "results"
params.min_length = 36
params.quality_threshold = 20

// ============================================================================
// PROCESS 1: Initial Quality Check
// ============================================================================

process FASTQC_RAW {
    tag "FastQC on ${sample_id}"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*.{html,zip}"
    tuple val(sample_id), path("*.html"), emit: html
    
    script:
    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}

// ============================================================================
// PROCESS 2: Adapter Trimming
// ============================================================================

process TRIMMOMATIC {
    tag "Trimming ${sample_id}"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R*.fastq.gz")
    
    script:
    def (read1, read2) = reads
    """
    trimmomatic PE -threads ${task.cpus} \\
        ${read1} ${read2} \\
        ${sample_id}_trimmed_R1.fastq.gz ${sample_id}_unpaired_R1.fastq.gz \\
        ${sample_id}_trimmed_R2.fastq.gz ${sample_id}_unpaired_R2.fastq.gz \\
        ILLUMINACLIP:adapters.fa:2:30:10 \\
        LEADING:${params.quality_threshold} \\
        TRAILING:${params.quality_threshold} \\
        SLIDINGWINDOW:4:${params.quality_threshold} \\
        MINLEN:${params.min_length}
    """
}

// ============================================================================
// PROCESS 3: Post-Trim Quality Check
// ============================================================================

process FASTQC_TRIMMED {
    tag "FastQC on trimmed ${sample_id}"
    publishDir "${params.outdir}/fastqc_trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*.{html,zip}"
    tuple val(sample_id), path("*.html"), emit: html
    
    script:
    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}

// ============================================================================
// PROCESS 4: Calculate Read Statistics
// ============================================================================

process READ_STATS {
    tag "Stats for ${sample_id}"
    publishDir "${params.outdir}/stats", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "${sample_id}_stats.txt"
    
    script:
    """
    echo "Sample: ${sample_id}" > ${sample_id}_stats.txt
    echo "Files: ${reads}" >> ${sample_id}_stats.txt
    echo "Read count:" >> ${sample_id}_stats.txt
    zcat ${reads[0]} | echo \$((\$(wc -l)/4)) >> ${sample_id}_stats.txt
    echo "File sizes:" >> ${sample_id}_stats.txt
    ls -lh ${reads} >> ${sample_id}_stats.txt
    """
}

// ============================================================================
// PROCESS 5: Aggregate All QC Reports
// ============================================================================

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    input:
    path "*"
    
    output:
    path "multiqc_report.html"
    path "multiqc_data"
    
    script:
    """
    multiqc .
    """
}

// ============================================================================
// WORKFLOW: Orchestrate All Processes
// ============================================================================

workflow {
    // Print welcome message
    log.info """
    ╔═══════════════════════════════════════════════════════════╗
    ║  RNA-seq Quality Control Pipeline - Week 1 Project      ║
    ╚═══════════════════════════════════════════════════════════╝
    
    Input reads    : ${params.reads}
    Output dir     : ${params.outdir}
    Min length     : ${params.min_length}
    Quality thresh : ${params.quality_threshold}
    
    """.stripIndent()
    
    // ========================================================================
    // STEP 1: Create input channel from paired FASTQ files
    // ========================================================================
    
    reads_ch = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { sample_id, files -> 
            // Extract just the sample name (remove _R1/_R2 suffix)
            def clean_id = sample_id.replaceAll(/_R[12]$/, '')
            [clean_id, files]
        }
    
    // ========================================================================
    // STEP 2: Initial quality check on raw reads
    // ========================================================================
    
    fastqc_raw = FASTQC_RAW(reads_ch)
    
    // ========================================================================
    // STEP 3: Trim adapters and low-quality bases
    // ========================================================================
    
    trimmed_ch = TRIMMOMATIC(reads_ch)
    
    // ========================================================================
    // STEP 4: Quality check on trimmed reads
    // ========================================================================
    
    fastqc_trimmed = FASTQC_TRIMMED(trimmed_ch)
    
    // ========================================================================
    // STEP 5: Calculate read statistics
    // ========================================================================
    
    stats = READ_STATS(trimmed_ch)
    
    // ========================================================================
    // STEP 6: Aggregate all QC outputs with MultiQC
    // ========================================================================
    
    // Collect all QC outputs
    all_qc = fastqc_raw.html
        .map { id, html -> html }
        .mix(fastqc_trimmed.html.map { id, html -> html })
        .mix(stats)
        .collect()
    
    // Generate final report
    MULTIQC(all_qc)
    
    // ========================================================================
    // COMPLETION MESSAGE
    // ========================================================================
    
    workflow.onComplete {
        log.info """
        ╔═══════════════════════════════════════════════════════════╗
        ║  Pipeline Completed!                                     ║
        ╚═══════════════════════════════════════════════════════════╝
        
        Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Duration  : ${workflow.duration}
        Output    : ${params.outdir}
        
        """.stripIndent()
    }
}

// ============================================================================
// WORKFLOW METADATA
// ============================================================================

workflow.onError {
    log.error "Pipeline failed!"
    log.error "Error message: ${workflow.errorMessage}"
}
```

### Understanding the Integration

Let's break down how this uses Week 1 concepts:

**Day 2 - Groovy Syntax**:
```groovy
// String interpolation
"${sample_id}_trimmed_R1.fastq.gz"

// Map transformation
.map { sample_id, files -> [clean_id, files] }

// Tuple destructuring
def (read1, read2) = reads

// Maps for parameters
params.quality_threshold
```

**Day 3 - Processes**:
```groovy
// Every process has clear input/output
process FASTQC_RAW {
    input: tuple val(sample_id), path(reads)
    output: path "*.html"
    script: "fastqc ${reads}"
}

// Directives like tag, publishDir
tag "FastQC on ${sample_id}"
publishDir "${params.outdir}/fastqc_raw"
```

**Day 4 - Channels**:
```groovy
// Queue channel from paired files
reads_ch = Channel.fromFilePairs(params.reads)

// Channel operators
.map { ... }
.mix(...)
.collect()

// Named outputs
emit: html
```

**Day 5 - Workflows**:
```groovy
// Sequential data flow
fastqc_raw = FASTQC_RAW(reads_ch)
trimmed_ch = TRIMMOMATIC(reads_ch)
fastqc_trimmed = FASTQC_TRIMMED(trimmed_ch)

// Branching (reads_ch used by multiple processes)
// Merging (mix and collect)
```

**Day 6 - Execution Understanding**:
```groovy
// When this runs:
// - Each task in work/xx/yyyy/
// - Can use -resume
// - publishDir saves results
// - Logs to .nextflow.log
```

### How to "Run" This (Conceptually)

Since we're not installing Nextflow, here's what would happen:

```bash
# Execution command
nextflow run qc_pipeline.nf \
  --reads "data/*_R{1,2}.fastq.gz" \
  --outdir results \
  -resume \
  -with-report report.html

# What happens:
# 1. Nextflow parses the workflow
# 2. Creates channels from input files
# 3. Executes processes in parallel
# 4. Each task in isolated directory
# 5. Collects outputs to results/
# 6. Generates MultiQC report
```

**Expected Timeline** (for 10 samples with 4 cores):
```
Time 0:00 - FASTQC_RAW starts (4 samples parallel)
Time 0:05 - FASTQC_RAW continues, TRIMMOMATIC starts
Time 0:10 - FASTQC_TRIMMED starts, READ_STATS starts
Time 0:15 - All samples processed, MULTIQC starts
Time 0:16 - Pipeline complete!

Total: ~16 minutes vs. ~100 minutes sequential
```

### Project Enhancements

Want to extend this? Try adding:

**Enhancement 1**: Add error handling
```groovy
process TRIMMOMATIC {
    errorStrategy 'retry'
    maxRetries 3
    
    // ... rest
}
```

**Enhancement 2**: Add resource specifications
```groovy
process FASTQC_RAW {
    cpus 2
    memory '4 GB'
    time '1h'
    
    // ... rest
}
```

**Enhancement 3**: Add conditional execution
```groovy
process DEEP_ANALYSIS {
    when:
    params.deep_mode == true
    
    // ... rest
}
```

---

## 🤔 Self-Assessment and Reflection (4 minutes)

### Knowledge Check

Rate your confidence (1-5) for each concept:

**Conceptual Understanding**:
- [ ] I understand what Nextflow is and when to use it (Day 1)
- [ ] I can explain orchestration vs processing (Day 1)
- [ ] I understand automatic parallelization (Day 1)

**Groovy Syntax**:
- [ ] I can write string interpolation with `${}` (Day 2)
- [ ] I can use `.collect {}` and `.findAll {}` (Day 2)
- [ ] I understand closures and `it` (Day 2)

**Processes**:
- [ ] I can write a complete process (Day 3)
- [ ] I understand input/output types (Day 3)
- [ ] I know when to use tuples (Day 3)

**Channels**:
- [ ] I understand queue vs value channels (Day 4)
- [ ] I can create channels from files (Day 4)
- [ ] I know why channels enable parallelization (Day 4)

**Workflows**:
- [ ] I can connect processes into workflows (Day 5)
- [ ] I understand data flow (Day 5)
- [ ] I know when to use `.collect()` (Day 5)

**Execution**:
- [ ] I understand the work directory (Day 6)
- [ ] I know how resume works (Day 6)
- [ ] I can debug failed tasks (Day 6)

### Areas for Review

**If you rated anything 1-2**, review that day's material.

**Common challenge areas**:
1. **Queue vs Value channels** - Review Day 4 if unclear
2. **When to use `.collect()`** - Review Day 5 examples
3. **Groovy syntax** - Review Day 2 exercises
4. **Process inputs/outputs** - Review Day 3 patterns

### Your Learning Style

**Reflect on how you learn best**:
- Do you prefer reading explanations or seeing examples?
- Do exercises help cement concepts?
- Do comparisons to Python help or confuse?
- Do you need more or less detail?

This helps you approach Week 2 optimally!

### Week 1 Achievements

Take a moment to appreciate what you've accomplished:

✅ **Conceptual foundation** - You understand Nextflow's purpose  
✅ **Language skills** - You can read and write Groovy  
✅ **Building blocks** - You can create processes  
✅ **Data flow** - You understand channels  
✅ **Orchestration** - You can build workflows  
✅ **Execution model** - You know how it runs  

**You've built a solid foundation in just 7 days!** 🎉

---

## 🎯 Looking Ahead: Week 2 Preview (2 minutes)

Next week, you'll add **practical skills** to your foundation:

**Week 2: Practical Pipeline Construction**

**Day 8**: Parameters and Workflow Flexibility
- Make workflows configurable
- Use `params` effectively
- Create reusable pipelines

**Day 9**: Working with Collections
- Advanced channel operators
- Data manipulation patterns
- Complex transformations

**Day 10**: Combining Multiple Inputs
- Handling paired data
- Using `each` qualifier
- Cross-product patterns

**Day 11**: Publishing and Managing Outputs
- `publishDir` directive
- Output organization
- Modes: copy, symlink, move

**Day 12**: Error Handling and Robustness
- Error strategies
- Retry mechanisms
- Input validation

**Day 13**: Working with Containers
- Docker/Singularity integration
- Reproducibility
- Per-process containers

**Day 14**: Debugging Workflows Systematically
- Advanced debugging techniques
- Week 2 integration project

**What makes Week 2 different**:
- Week 1: Building blocks and understanding
- Week 2: Making pipelines robust and production-ready
- Week 3: Advanced patterns and optimization
- Week 4: Complete production workflows

---

## ✅ Week 1 Completion Checklist

Before moving to Week 2, ensure you:

**Understanding**:
- [ ] Can explain Nextflow's value proposition
- [ ] Know when to use Nextflow vs Python
- [ ] Understand automatic parallelization

**Syntax**:
- [ ] Can write Groovy string interpolation
- [ ] Can use closures with `it`
- [ ] Can work with Lists and Maps

**Processes**:
- [ ] Can write processes with input/output/script
- [ ] Understand process isolation
- [ ] Know all input/output types

**Channels**:
- [ ] Understand queue vs value channels
- [ ] Can create channels from various sources
- [ ] Know how channels enable parallelization

**Workflows**:
- [ ] Can connect processes into workflows
- [ ] Understand data flow patterns
- [ ] Can use channel operators

**Execution**:
- [ ] Understand work directory structure
- [ ] Know how resume works
- [ ] Can debug failed tasks

**Integration**:
- [ ] Can read and understand the Week 1 project
- [ ] Could modify it for different tools
- [ ] Understand how all pieces fit together

---

## 🎊 Congratulations!

**You've completed Week 1 of your Nextflow journey!**

In just 7 days, you've gone from zero Nextflow knowledge to understanding how to build complete bioinformatics workflows. That's a significant achievement!

### What You Can Do Now

✅ Read and understand Nextflow workflows  
✅ Write basic processes  
✅ Work with channels  
✅ Build simple pipelines  
✅ Understand execution  
✅ Debug issues  

### What's Next

Week 2 will teach you to make your workflows:
- **Flexible** with parameters
- **Robust** with error handling
- **Reproducible** with containers
- **Professional** with proper output management

### Keep Going!

The hardest part is behind you - **understanding the fundamentals**. Everything from here builds on what you know.

Take a break, celebrate your progress, and come back ready for Week 2! 🚀

---

## 📝 Week 1 Summary Card

**Save this for quick reference**:

```
WEEK 1: FOUNDATIONS
═══════════════════

DAY 1: Orchestration vs Processing
  → Nextflow coordinates, Python processes

DAY 2: Groovy Essentials
  → "${var}", .collect{}, .findAll{}

DAY 3: Processes
  → input/output/script blocks
  → val, path, tuple

DAY 4: Channels
  → Queue (once), Value (reusable)
  → fromPath(), value()

DAY 5: Workflows
  → process outputs → channel inputs
  → .collect() for aggregation

DAY 6: Execution
  → work/xx/yyyy/ per task
  → -resume uses task hashing

INTEGRATION: Complete QC Pipeline
  → All concepts in one workflow
```

---

**Completed Week 1?** Update your `PROGRESS.md` with your achievements! 🎉

**Your progress**: 7/28 days (25%) complete - **Quarter done!**

**Tomorrow**: Week 2 begins! Rest up and prepare for practical pipeline skills! 🚀

---

*Week 2 Day 8: Parameters and Workflow Flexibility awaits!*

**See you in Week 2! You've got this! 🎊**