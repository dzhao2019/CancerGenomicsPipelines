# Day 5: Subworkflows and Production-Ready Pipelines

**Learning Time**: 30 minutes  
**Prerequisites**: Completed Days 1-4  
**Goal**: Organize workflows for production use and build reusable components

---

## 📖 Introduction (2 minutes)

Welcome to Day 5! You're on the verge of professional pipeline development.

Over the last 4 days, you've learned:
- **Day 1:** What Nextflow is
- **Day 2:** Groovy syntax
- **Day 3:** How to write processes
- **Day 4:** How to connect them with workflows

Today you'll learn how to **organize workflows professionally** so they're:
- Reusable across projects
- Easy to maintain
- Simple to modify
- Ready for production use

### The Challenge You're Solving

Imagine you have a 10-process pipeline. Each process definition is long. Your main workflow is complex. It's hard to find things. It's hard to reuse parts. It's messy.

**Solution: Subworkflows**

A subworkflow is a workflow within a workflow. It groups related processes together, making your code:
- Modular (small, reusable pieces)
- Readable (easier to understand)
- Maintainable (easier to update)
- Professional (organized like production code)

---

## 🎯 Learning Objectives

By the end of Day 5, you should be able to:

1. **Explain subworkflows** - What they are and why they matter
2. **Create subworkflows** - Group related processes
3. **Reuse subworkflows** - Use the same subworkflow in different ways
4. **Organize code** - Structure for production pipelines
5. **Configure pipelines** - Use parameters and config files
6. **Understand best practices** - How professionals organize Nextflow code
7. **Build modular pipelines** - Combine subworkflows into larger workflows

---

## 📚 Key Concepts (20 minutes)

### 1. What is a Subworkflow?

A **subworkflow** is a workflow that's called from another workflow.

**Simple example:**

```groovy
// Regular workflow (what you've been writing)
workflow {
    samples = Channel.fromPath("*.fastq")
    PROCESS1(samples)
    PROCESS2(samples)
}

// Subworkflow (a workflow that can be called)
workflow quality_control {
    take: samples  // Input
    
    main:
    FASTQC(samples)
    TRIMMOMATIC(samples)
    
    emit: // Output
    TRIMMOMATIC.out
}

// Main workflow that uses the subworkflow
workflow {
    samples = Channel.fromPath("*.fastq")
    qc_results = quality_control(samples)  // Call subworkflow!
}
```

**Key differences:**
- Subworkflow has `take` block (what it receives)
- Subworkflow has `emit` block (what it outputs)
- Main workflow can call it like a process

### Why Use Subworkflows?

#### Problem 1: Long, Complex Code
```groovy
// Without subworkflows - hard to read
workflow {
    samples = Channel.fromPath("*.fastq")
    
    // QC section (5 processes)
    qc_raw = FASTQC(samples)
    trimmed = TRIMMOMATIC(samples)
    qc_trimmed = FASTQC_POST(trimmed)
    
    // Alignment section (3 processes)
    aligned = ALIGN(trimmed)
    sorted = SORT_BAM(aligned)
    indexed = INDEX_BAM(sorted)
    
    // Variant calling section (2 processes)
    variants = CALL_VARIANTS(indexed)
    filtered = FILTER_VCF(variants)
    
    filtered.view()
}
```

**Solution: Subworkflows**
```groovy
// With subworkflows - clear organization
workflow {
    samples = Channel.fromPath("*.fastq")
    
    qc_out = quality_control(samples)
    alignment_out = alignment(qc_out)
    variants_out = variant_calling(alignment_out)
    
    variants_out.view()
}
```

#### Problem 2: Code Reuse
Without subworkflows, if you want to use QC in another project, you copy-paste processes.

With subworkflows, you reuse the entire `quality_control` subworkflow.

#### Problem 3: Testing
It's hard to test individual pieces of a monolithic workflow.

With subworkflows, you can test each piece separately.

---

### 2. Anatomy of a Subworkflow

```groovy
workflow quality_control {
    // 1. TAKE block: what this subworkflow receives
    take:
        fastq_files
    
    // 2. MAIN block: the work
    main:
        raw_qc = FASTQC(fastq_files)
        trimmed = TRIMMOMATIC(fastq_files)
        trimmed_qc = FASTQC_POST(trimmed)
    
    // 3. EMIT block: what this subworkflow outputs
    emit:
        reports: raw_qc          // Can name outputs
        trimmed_reads: trimmed
        trim_reports: trimmed_qc
}
```

**How to call it:**
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    
    // Call subworkflow
    qc = quality_control(samples)
    
    // Access named outputs
    qc.reports.view()
    qc.trimmed_reads.view()
}
```

---

### 3. Real Example: Multi-Stage Pipeline

Let's build a realistic pipeline with three subworkflows:

```groovy
// Subworkflow 1: Quality Control
workflow quality_control {
    take:
        reads
    
    main:
        fastqc_raw = FASTQC(reads)
        trimmed = TRIMMOMATIC(reads)
        fastqc_trim = FASTQC_POST(trimmed)
    
    emit:
        trimmed_reads: trimmed
        qc_reports: fastqc_trim.mix(fastqc_raw)
}

// Subworkflow 2: Alignment
workflow alignment {
    take:
        reads
        reference
    
    main:
        aligned = ALIGN(reads, reference)
        sorted = SORT_BAM(aligned)
        indexed = INDEX_BAM(sorted)
    
    emit:
        bam_files: indexed
}

// Subworkflow 3: Variant Calling
workflow variant_calling {
    take:
        bam_files
        reference
    
    main:
        vcf = CALL_VARIANTS(bam_files, reference)
        filtered = FILTER_VCF(vcf)
    
    emit:
        variants: filtered
}

// Main workflow: orchestrates everything
workflow {
    // Input
    reads = Channel.fromPath("data/*.fastq")
        .map { f -> [f.baseName, f] }
    reference = Channel.value("reference/hg38.fa")
    
    // Quality control
    qc_out = quality_control(reads)
    
    // Alignment (use trimmed reads from QC)
    align_out = alignment(qc_out.trimmed_reads, reference)
    
    // Variant calling
    variants = variant_calling(align_out.bam_files, reference)
    
    // Output
    variants.view()
}
```

**Benefits:**
- Clear stages (QC → Alignment → Variants)
- Each subworkflow is independent
- Easy to modify one stage without affecting others
- Could reuse QC subworkflow in another pipeline

---

### 4. Parameters and Configuration

A professional pipeline is configurable:

```groovy
// params.nf or at top of main.nf
params {
    input_dir = "data"
    output_dir = "results"
    reference = "reference/hg38.fa"
    
    // Tool-specific parameters
    quality_threshold = 30
    min_depth = 10
    
    // Resource parameters
    alignment_cpus = 8
    alignment_memory = "16GB"
    
    // Feature flags
    skip_qc = false
    skip_trim = false
}
```

**Use in processes:**
```groovy
process ALIGN {
    cpus params.alignment_cpus
    memory params.alignment_memory
    
    input:
        tuple val(id), path(fastq)
        path reference
    output:
        tuple val(id), path("${id}.bam")
    script:
        """
        bwa mem -t ${task.cpus} ${reference} ${fastq} > ${id}.bam
        """
}
```

**Use in workflow:**
```groovy
workflow {
    samples = Channel.fromPath("${params.input_dir}/*.fastq")
    reference = Channel.value(params.reference)
    
    // Use parameters to control flow
    if (!params.skip_qc) {
        qc_out = quality_control(samples)
    }
    
    align_out = alignment(samples, reference)
}
```

**Command line override:**
```bash
nextflow run main.nf \
    --input_dir data/new \
    --reference reference/hg19.fa \
    --quality_threshold 20
```

---

### 5. File Organization: The Professional Structure

For a production pipeline, organize files like this:

```
my_pipeline/
├── main.nf                 # Main workflow
├── nextflow.config         # Configuration
├── params.json             # Default parameters
├── modules/
│   ├── qc.nf              # QC subworkflow + processes
│   ├── alignment.nf       # Alignment subworkflow + processes
│   ├── variants.nf        # Variant calling subworkflow + processes
│   └── utils.nf           # Helper functions
├── workflows/
│   ├── qc_only.nf         # Just QC workflow
│   └── full_pipeline.nf   # Complete pipeline
├── conf/
│   ├── base.config        # Default resources
│   ├── local.config       # Local machine settings
│   ├── cluster.config     # HPC cluster settings
│   └── docker.config      # Container settings
├── tests/
│   ├── test_data/
│   └── test_modules.nf    # Module unit tests
└── docs/
    ├── README.md
    └── USAGE.md
```

### Importing Modules

Instead of everything in one file, use `include`:

```groovy
// main.nf
include { FASTQC; TRIMMOMATIC } from './modules/qc.nf'
include { ALIGN; SORT_BAM } from './modules/alignment.nf'

workflow {
    samples = Channel.fromPath("*.fastq")
    FASTQC(samples)
    // Processes are now available!
}
```

### Using nf-core (Industry Standard)

The professional way uses `nf-core`:

```groovy
// Include a subworkflow from nf-core
include { FASTQC_CHECK } from './subworkflows/nf-core/fastq_fastqc_umitools_fastp/main'

workflow {
    // Uses pre-built, tested subworkflows
    FASTQC_CHECK(samples)
}
```

---

### 6. Conditional Execution: Skip Steps

```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    reference = Channel.value("ref.fa")
    
    // Skip QC if requested
    if (params.skip_qc) {
        qc_samples = samples
    } else {
        qc_out = quality_control(samples)
        qc_samples = qc_out.trimmed
    }
    
    // Always do alignment
    alignment(qc_samples, reference)
}
```

---

### 7. Combining Subworkflows: Building Larger Workflows

You can have subworkflows that call other subworkflows:

```groovy
// Small subworkflows
workflow preprocessing {
    take: reads
    main:
        qc_out = quality_control(reads)
    emit:
        trimmed: qc_out.trimmed
}

workflow processing {
    take: reads
    main:
        align_out = alignment(reads)
    emit:
        bam: align_out.bam
}

// Larger subworkflow that combines them
workflow analysis {
    take: reads
    main:
        preprocess_out = preprocessing(reads)
        process_out = processing(preprocess_out.trimmed)
    emit:
        results: process_out.bam
}

// Main workflow
workflow {
    samples = Channel.fromPath("*.fastq")
    analysis(samples)
}
```

---

### 8. Error Handling in Production

```groovy
workflow {
    samples = Channel.fromPath("data/*.fastq")
        .ifEmpty { error("No FASTQ files found!") }
    
    reference = Channel.value("reference.fa")
        .ifEmpty { error("Reference not found!") }
    
    // Check parameters
    if (!params.output_dir) {
        error("output_dir parameter required!")
    }
    
    // Run workflow with error handling
    try {
        qc_out = quality_control(samples)
        align_out = alignment(qc_out, reference)
    } catch (Exception e) {
        log.error "Pipeline failed: ${e.message}"
        System.exit(1)
    }
}
```

---

### 9. Logging and Progress Reports

```groovy
workflow {
    log.info """
    ╔═══════════════════════════════════════════════════════════╗
    ║         RNA-Seq Analysis Pipeline                         ║
    ╠═══════════════════════════════════════════════════════════╣
    ║ Input Directory : ${params.input_dir}
    ║ Reference      : ${params.reference}
    ║ Output Dir     : ${params.output_dir}
    ║ Skip QC        : ${params.skip_qc}
    ╚═══════════════════════════════════════════════════════════╝
    """.stripIndent()
    
    samples = Channel.fromPath("${params.input_dir}/*.fastq")
    
    log.info "Found ${samples.size()} samples"
    
    qc_out = quality_control(samples)
    
    log.info "Quality control completed"
}
```

---

### 10. Best Practices Summary

**DO:**
✅ Use subworkflows to group related processes  
✅ Use parameters for configuration  
✅ Organize code into modules  
✅ Add error handling  
✅ Include logging  
✅ Document your pipeline  
✅ Use consistent naming (PROCESS_NAME, workflow_name)  
✅ Keep processes focused (one job each)  
✅ Preserve metadata (sample IDs) through pipeline  

**DON'T:**
❌ Put everything in one file  
❌ Hardcode paths and parameters  
❌ Ignore errors  
❌ Skip logging  
❌ Create overly complex workflows  
❌ Mix concerns (QC and alignment in one process)  

---

## 🔗 Connecting to Previous Days

### Days 1-3: Foundation
- Day 1: What Nextflow is
- Day 2: Groovy syntax
- Day 3: How to write processes

### Day 4: Workflows
- How to connect processes

### Day 5: Organization
- How to organize workflows professionally

**Together:** You can build production-ready pipelines!

---

## ✅ Completion Checklist

- [ ] I understand what a subworkflow is
- [ ] I can create a subworkflow
- [ ] I can use take/main/emit blocks
- [ ] I know when to use subworkflows
- [ ] I understand parameters and configuration
- [ ] I know how to organize code in modules
- [ ] I can import processes from other files
- [ ] I understand best practices
- [ ] I can build modular, professional pipelines

---

## 🔑 Key Takeaways

### Subworkflows
A subworkflow is a workflow within a workflow that groups related processes together.

```groovy
workflow my_subworkflow {
    take: input_data
    main:
        PROCESS1(input_data)
    emit:
        PROCESS1.out
}
```

### Parameters
Make pipelines configurable with parameters instead of hardcoding values.

```groovy
params {
    input = "data/"
    reference = "ref.fa"
}
```

### Modular Organization
Split code across files for easier maintenance and reuse.

```groovy
include { PROCESS1 } from './modules/module1.nf'
```

### Professional Structure
Organize files logically with modules, workflows, config, tests, and docs.

---

## 🚀 Ready for Exercises?

You now understand:
- Subworkflows ✅
- Parameters ✅
- Modular organization ✅
- Best practices ✅

Time to build production-ready pipelines!

---

*This is Day 5 of 28. You're now ready to build professional bioinformatics pipelines!*
