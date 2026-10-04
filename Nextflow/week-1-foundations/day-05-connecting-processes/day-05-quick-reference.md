# Day 5: Subworkflows and Production Pipelines - Quick Reference

**Print this for quick lookup while organizing workflows!**

---

## 🏗️ Subworkflow Template

```groovy
workflow subworkflow_name {
    take:
        input_channel
    
    main:
        PROCESS1(input_channel)
        result = PROCESS2(input_channel)
    
    emit:
        output_name: result
}

// Use it
workflow {
    data = Channel.fromPath("*.fastq")
    out = subworkflow_name(data)
    out.output_name.view()
}
```

---

## 📋 Subworkflow Structure (3 Blocks)

### Block 1: TAKE (Input)
```groovy
take:
    channel1          // Input channel
    channel2          // Another input
```
What the subworkflow receives.

### Block 2: MAIN (Work)
```groovy
main:
    result = PROCESS(channel1, channel2)
    final = PROCESS2(result)
```
The actual processes and workflow logic.

### Block 3: EMIT (Output)
```groovy
emit:
    results: result        // Named output
    final_out: final       // Another named output
```
What the subworkflow outputs.

---

## 🔄 Using Subworkflows

### Simple Call
```groovy
workflow {
    data = Channel.fromPath("*.fastq")
    output = my_subworkflow(data)
}
```

### Multiple Inputs
```groovy
workflow {
    reads = Channel.fromPath("*.fastq")
    reference = Channel.value("ref.fa")
    
    output = my_subworkflow(reads, reference)
}
```

### Access Named Outputs
```groovy
workflow {
    data = Channel.fromPath("*.fastq")
    qc = quality_control(data)
    
    qc.reports.view()      // Access named output
    qc.trimmed.view()      // Another named output
}
```

### Chain Subworkflows
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    
    prep = preprocessing(samples)
    align = alignment(prep.trimmed_reads)
    results = variant_calling(align.bam_files)
}
```

---

## 📋 Parameters (Configuration)

### Define Parameters
```groovy
params {
    input_dir = "data"
    output_dir = "results"
    reference = "reference/hg38.fa"
    threads = 8
    skip_qc = false
}
```

### Use in Process
```groovy
process ALIGN {
    cpus params.threads
    memory "${params.memory}GB"
    
    script:
        """
        bwa mem -t ${task.cpus} ...
        """
}
```

### Use in Workflow
```groovy
workflow {
    samples = Channel.fromPath("${params.input_dir}/*.fastq")
    
    if (!params.skip_qc) {
        qc_out = quality_control(samples)
    }
}
```

### Override from Command Line
```bash
nextflow run main.nf \
    --input_dir data/new \
    --threads 16 \
    --skip_qc true
```

---

## 🗂️ File Organization

### Minimal Structure
```
my_pipeline/
├── main.nf
└── modules/
    └── processes.nf
```

### Standard Structure
```
my_pipeline/
├── main.nf
├── nextflow.config
├── modules/
│   ├── preprocessing.nf
│   ├── alignment.nf
│   └── variants.nf
└── workflows/
    ├── full.nf
    └── qc_only.nf
```

### Professional Structure
```
my_pipeline/
├── main.nf
├── nextflow.config
├── params.json
├── modules/
│   ├── preprocessing.nf
│   ├── alignment.nf
│   └── variants.nf
├── workflows/
│   └── main.nf
├── conf/
│   ├── base.config
│   ├── local.config
│   └── cluster.config
├── tests/
│   └── test_main.nf
└── docs/
    ├── README.md
    └── USAGE.md
```

---

## 📦 Importing Code

### Include Processes
```groovy
include { PROCESS1; PROCESS2 } from './modules/module.nf'
```

### Include Workflows
```groovy
include { my_workflow } from './modules/subworkflows.nf'
```

### Include with Alias
```groovy
include { PROCESS1 as ALIGN_READS } from './modules/alignment.nf'

// Use as
ALIGN_READS(reads)
```

### Import Multiple from One File
```groovy
include { 
    FASTQC
    TRIMMOMATIC
    FASTQC_POST
} from './modules/preprocessing.nf'
```

---

## ✅ Workflow Checklist

- [ ] Parameters defined with defaults
- [ ] Input validation (files exist?)
- [ ] Error handling (try/catch)
- [ ] Logging at key points
- [ ] Help/usage message
- [ ] Named outputs in subworkflows
- [ ] Code organized into modules
- [ ] Configuration file for settings

---

## 🔍 Error Handling

```groovy
// Check file exists
if (!file(params.reference).exists()) {
    error("Reference not found: ${params.reference}")
}

// Check channel not empty
Channel.fromPath("*.fastq")
    .ifEmpty { error("No FASTQ files found") }

// Try/catch for workflow
try {
    qc = quality_control(samples)
    align = alignment(qc.out)
} catch (Exception e) {
    log.error "Pipeline failed: ${e.message}"
    System.exit(1)
}

// Conditional parameter check
if (!params.output_dir) {
    error("output_dir parameter required")
}
```

---

## 📝 Logging and Messages

```groovy
// Info messages
log.info "Starting pipeline..."
log.info "Found ${samples.size()} samples"

// Warning messages
log.warn "Low quality threshold: ${params.quality}"

// Error messages
log.error "Critical error occurred"

// Beautiful header
log.info """
╔═══════════════════════════════════════════════╗
║              RNA-Seq Pipeline v1.0           ║
╠═══════════════════════════════════════════════╣
║ Input: ${params.input_dir}
║ Output: ${params.output_dir}
╚═══════════════════════════════════════════════╝
""".stripIndent()
```

---

## 🏗️ Common Patterns

### Pattern: Skip Optional Steps
```groovy
if (params.skip_qc) {
    processed = input_reads
} else {
    qc_out = quality_control(input_reads)
    processed = qc_out.trimmed
}

aligned = alignment(processed)
```

### Pattern: Multiple Subworkflows
```groovy
prep = preprocessing(samples)
align = alignment(prep.out)
variants = variant_calling(align.out)
counts = counting(align.bam)
```

### Pattern: Named Outputs
```groovy
workflow quality_control {
    take: reads
    main:
        raw = FASTQC(reads)
        trimmed = TRIM(reads)
        post = FASTQC(trimmed)
    emit:
        raw_qc: raw
        trimmed_reads: trimmed
        trimmed_qc: post
}
```

### Pattern: Conditional Output
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    
    qc = quality_control(samples)
    
    // Different next step based on parameter
    if (params.alignment_tool == "bwa") {
        aligned = bwa_align(qc.trimmed)
    } else if (params.alignment_tool == "bowtie2") {
        aligned = bowtie2_align(qc.trimmed)
    }
}
```

---

## 📊 Example: Complete Subworkflow

```groovy
workflow preprocessing {
    take:
        fastq_files
    
    main:
        // QC on raw reads
        raw_qc = FASTQC(fastq_files)
        
        // Trim adapters
        trimmed = TRIMMOMATIC(fastq_files)
        
        // QC on trimmed reads
        trim_qc = FASTQC_POST(trimmed)
    
    emit:
        trimmed_reads: trimmed
        qc_reports: raw_qc.mix(trim_qc)
}

// Use it
workflow {
    reads = Channel.fromPath("data/*.fastq")
        .map { f -> [f.baseName, f] }
    
    preprocess = preprocessing(reads)
    
    // Trimmed reads for next step
    aligned = ALIGN(preprocess.trimmed_reads)
    
    aligned.view()
}
```

---

## 🎯 Best Practices

**DO:**
✅ Use subworkflows to group related processes  
✅ Use parameters instead of hardcoding  
✅ Add error handling and validation  
✅ Include logging for debugging  
✅ Organize code into modules  
✅ Use clear, descriptive names  
✅ Document your parameters  
✅ Version your pipeline  

**DON'T:**
❌ Put everything in one file  
❌ Hardcode paths or settings  
❌ Ignore errors  
❌ Skip logging  
❌ Create overly complex subworkflows  
❌ Use vague names (workflow, process1)  
❌ Copy-paste code between projects  

---

## 📋 Quick Start Templates

### Simple Pipeline
```groovy
params {
    input = "data/*.fastq"
    output = "results"
}

workflow {
    samples = Channel.fromPath(params.input)
    PROCESS(samples)
}
```

### Two-Step Pipeline
```groovy
workflow step1 {
    take: input
    main:
        PROC1(input)
    emit:
        PROC1.out
}

workflow step2 {
    take: input
    main:
        PROC2(input)
    emit:
        PROC2.out
}

workflow {
    samples = Channel.fromPath("*.fastq")
    s1 = step1(samples)
    s2 = step2(s1)
}
```

### With Error Handling
```groovy
params {
    input = "data"
    output = "results"
}

workflow {
    input = Channel.fromPath("${params.input}/*.fastq")
        .ifEmpty { error("No input files") }
    
    log.info "Processing ${input.size()} samples"
    
    try {
        PROCESS(input)
    } catch (e) {
        log.error "Failed: ${e}"
        System.exit(1)
    }
}
```

---

## 🔗 Configuration Files

### nextflow.config
```groovy
// Set defaults
params {
    input = "data"
    output = "results"
}

// Process settings
process {
    publishDir = params.output
}

// Profiles
profiles {
    local {
        executor = 'local'
    }
    cluster {
        executor = 'slurm'
        queue = 'normal'
    }
}
```

### Using Config
```bash
nextflow run main.nf -c nextflow.config
nextflow run main.nf -profile cluster
```

---

*Print this page and keep it handy while building production pipelines!*

**Next:** Day 6 - Advanced patterns and practical techniques

