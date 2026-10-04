# Day 5: Subworkflows and Production Pipelines - Hands-On Exercises

**Time**: 20 minutes for exercises + 15 minutes for solutions  
**Difficulty**: Moderate to Advanced  
**Goal**: Build and organize production-ready workflows

---

## Exercise 1: Understanding Subworkflows (5 minutes)

### Part A: Identify Subworkflow Components

**Given this subworkflow, label each part:**

```groovy
workflow quality_control {
    take:
        fastq_files
    
    main:
        raw_qc = FASTQC(fastq_files)
        trimmed = TRIMMOMATIC(fastq_files)
        trim_qc = FASTQC_POST(trimmed)
    
    emit:
        trimmed_reads: trimmed
        reports: trim_qc
}
```

**Identify:**
1. The `take` block and what it receives
2. The `main` block and what it does
3. The `emit` block and what it outputs
4. How many outputs does this subworkflow have?

<details>
<summary>✓ Solution A</summary>

**1. The take block:**
```groovy
take:
    fastq_files
```
- Receives: FASTQ files channel
- This is the input to the subworkflow

**2. The main block:**
```groovy
main:
    raw_qc = FASTQC(fastq_files)
    trimmed = TRIMMOMATIC(fastq_files)
    trim_qc = FASTQC_POST(trimmed)
```
- Runs FASTQC on raw reads
- Trims adapters with TRIMMOMATIC
- Runs FASTQC again on trimmed reads
- This is the work of the subworkflow

**3. The emit block:**
```groovy
emit:
    trimmed_reads: trimmed
    reports: trim_qc
```
- Outputs the trimmed reads
- Outputs the QC reports
- These become available to calling workflow

**4. How many outputs?**
Two named outputs:
- `trimmed_reads` - the processed FASTQ files
- `reports` - the QC HTML reports

**How to use it:**
```groovy
workflow {
    reads = Channel.fromPath("*.fastq")
    qc = quality_control(reads)
    
    // Access outputs by name
    qc.trimmed_reads.view()  // Trimmed FASTQ
    qc.reports.view()         // QC reports
}
```

</details>

---

### Part B: Compare: Process vs Subworkflow

**How is a subworkflow different from a process?**

| Aspect | Process | Subworkflow |
|--------|---------|-------------|
| What it contains | ??? | ??? |
| Input declaration | `input:` | ??? |
| Output declaration | `output:` | ??? |
| Internal work | Bash/Python script | ??? |
| Used for | ??? | ??? |

<details>
<summary>✓ Solution B</summary>

| Aspect | Process | Subworkflow |
|--------|---------|-------------|
| What it contains | Single tool/command | Multiple processes |
| Input declaration | `input:` | `take:` |
| Output declaration | `output:` | `emit:` |
| Internal work | Bash/Python script | Groovy workflow logic |
| Used for | Running a single tool | Grouping related processes |

**Key difference:**
- **Process:** Wraps one tool/command
- **Subworkflow:** Orchestrates multiple processes

**Analogy:**
- Process = a function that does one thing
- Subworkflow = a function that calls other functions

</details>

---

## Exercise 2: Write a Simple Subworkflow (10 minutes)

### Your Task

**Convert this workflow into a subworkflow:**

```groovy
// Current workflow
workflow {
    samples = Channel.fromPath("data/*.fastq")
    
    qc1 = FASTQC(samples)
    trimmed = TRIMMOMATIC(samples)
    qc2 = FASTQC_POST(trimmed)
}
```

**Requirements:**
- Create subworkflow named `preprocessing`
- Takes FASTQ files as input
- Outputs trimmed reads
- Call it from main workflow

**Your subworkflow:**

```groovy
workflow preprocessing {
    ???
}

workflow {
    ???
}
```

<details>
<summary>✓ Solution</summary>

```groovy
workflow preprocessing {
    take:
        fastq_files
    
    main:
        qc1 = FASTQC(fastq_files)
        trimmed = TRIMMOMATIC(fastq_files)
        qc2 = FASTQC_POST(trimmed)
    
    emit:
        trimmed_reads: trimmed
        qc_reports: qc1.mix(qc2)  // Combine both QC outputs
}

workflow {
    samples = Channel.fromPath("data/*.fastq")
    
    preprocess = preprocessing(samples)
    
    preprocess.trimmed_reads.view()
    preprocess.qc_reports.view()
}
```

**Key points:**
- `take:` declares input (`fastq_files`)
- `main:` contains the workflow logic
- `emit:` declares outputs with names
- `preprocess(samples)` calls the subworkflow
- Access outputs with dot notation: `preprocess.trimmed_reads`

**More complete version with QC outputs:**

```groovy
workflow preprocessing {
    take:
        fastq_files
    
    main:
        qc_raw = FASTQC(fastq_files)
        trimmed = TRIMMOMATIC(fastq_files)
        qc_trimmed = FASTQC_POST(trimmed)
    
    emit:
        trimmed_reads: trimmed
        raw_reports: qc_raw
        trimmed_reports: qc_trimmed
}

workflow {
    samples = Channel.fromPath("data/*.fastq")
    
    preprocess = preprocessing(samples)
    
    // Use trimmed reads for next step
    aligned = ALIGN(preprocess.trimmed_reads)
    
    aligned.view()
}
```

</details>

---

## Exercise 3: Multi-Level Subworkflows (10 minutes)

### Your Task

**Build a pipeline with multiple subworkflows:**

You have three groups of processes:
1. **Preprocessing:** FASTQC → TRIMMOMATIC → FASTQC_POST
2. **Alignment:** ALIGN → SORT_BAM → INDEX_BAM
3. **Variant Calling:** CALL_VARIANTS → FILTER_VCF

Create three subworkflows and a main workflow that chains them together.

**Your code:**

```groovy
workflow preprocessing {
    ???
}

workflow alignment {
    ???
}

workflow variant_calling {
    ???
}

workflow {
    ???
}
```

<details>
<summary>✓ Solution</summary>

```groovy
workflow preprocessing {
    take:
        reads
    
    main:
        FASTQC(reads)
        trimmed = TRIMMOMATIC(reads)
        FASTQC_POST(trimmed)
    
    emit:
        trimmed_reads: trimmed
}

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

workflow {
    // Input
    reads = Channel.fromPath("data/*.fastq")
        .map { f -> [f.baseName, f] }
    reference = Channel.value("reference/hg38.fa")
    
    // Chain subworkflows
    prep = preprocessing(reads)
    align = alignment(prep.trimmed_reads, reference)
    variants = variant_calling(align.bam_files, reference)
    
    // Output
    variants.view()
}
```

**Data flow:**
```
reads → preprocessing → trimmed_reads
                          ↓
                       alignment → bam_files
                                      ↓
                             variant_calling → variants
```

**Key insight:** Each subworkflow is independent. You can:
- Reuse `preprocessing` in other pipelines
- Skip a stage if needed
- Test each stage separately

</details>

---

## Exercise 4: Using Parameters (10 minutes)

### Part A: Add Parameters

**Add configuration to this workflow:**

```groovy
workflow {
    samples = Channel.fromPath("data/*.fastq")
    reference = Channel.value("reference/hg38.fa")
    
    qc = quality_control(samples)
    aligned = alignment(qc.trimmed, reference)
}
```

**Requirements:**
- Make `data` directory configurable
- Make reference path configurable
- Make alignment CPU count configurable
- Set defaults

<details>
<summary>✓ Solution A</summary>

```groovy
// Add at top of file
params {
    input_dir = "data"
    reference = "reference/hg38.fa"
    align_cpus = 8
}

// Use in processes
process ALIGN {
    cpus params.align_cpus
    
    input:
        tuple val(id), path(fastq)
        path ref
    output:
        tuple val(id), path("${id}.bam")
    script:
        """
        bwa mem -t ${task.cpus} ${ref} ${fastq} > ${id}.bam
        """
}

// Use in workflow
workflow {
    samples = Channel.fromPath("${params.input_dir}/*.fastq")
        .map { f -> [f.baseName, f] }
    reference = Channel.value(params.reference)
    
    qc = quality_control(samples)
    aligned = alignment(qc.trimmed, reference)
}
```

**Now you can override from command line:**
```bash
nextflow run main.nf \
    --input_dir data/new_samples \
    --reference reference/hg19.fa \
    --align_cpus 16
```

</details>

---

### Part B: Conditional Execution

**Add ability to skip QC:**

```groovy
params {
    input_dir = "data"
    skip_qc = false
}

workflow {
    samples = Channel.fromPath("${params.input_dir}/*.fastq")
    
    // Skip QC if requested
    if (params.skip_qc) {
        processed = samples
    } else {
        qc_out = quality_control(samples)
        processed = qc_out.trimmed
    }
    
    aligned = alignment(processed)
}
```

**Run without QC:**
```bash
nextflow run main.nf --skip_qc true
```

---

## Exercise 5: Organizing Code into Modules (10 minutes)

### Your Task

**Split code into multiple files:**

Current structure:
```
pipeline.nf    (all processes and workflows)
```

Desired structure:
```
main.nf              (main workflow)
modules/
├── preprocessing.nf (FASTQC, TRIMMOMATIC)
├── alignment.nf     (ALIGN, SORT_BAM, INDEX)
└── variants.nf      (CALL_VARIANTS, FILTER)
```

**In main.nf:**

```groovy
include { FASTQC; TRIMMOMATIC } from './modules/preprocessing'
include { ALIGN; SORT_BAM } from './modules/alignment'
include { CALL_VARIANTS } from './modules/variants'

workflow {
    samples = Channel.fromPath("*.fastq")
    
    FASTQC(samples)
    // Processes are now available!
}
```

**In modules/preprocessing.nf:**

```groovy
process FASTQC {
    input: tuple val(id), path(fastq)
    output: tuple val(id), path("*_fastqc.html")
    script: """
        fastqc ${fastq}
    """
}

process TRIMMOMATIC {
    input: tuple val(id), path(fastq)
    output: tuple val(id), path("*.trimmed.fastq")
    script: """
        trimmomatic SE ${fastq} output.trimmed.fastq ...
    """
}

// Can also include workflow
workflow preprocessing {
    take: reads
    main:
        FASTQC(reads)
        trimmed = TRIMMOMATIC(reads)
    emit:
        trimmed_reads: trimmed
}
```

**Key benefits:**
- Easier to find code
- Can reuse modules in other projects
- Easier to test individual modules
- Professional organization

---

## Exercise 6: Error Handling and Logging (10 minutes)

### Add Professional Error Handling

```groovy
params {
    input_dir = "data"
    reference = "reference/hg38.fa"
    output_dir = "results"
}

workflow {
    // Check inputs exist
    input = Channel.fromPath("${params.input_dir}/*.fastq")
        .ifEmpty { error("No FASTQ files found in ${params.input_dir}") }
    
    if (!file(params.reference).exists()) {
        error("Reference file not found: ${params.reference}")
    }
    
    if (!params.output_dir) {
        error("output_dir parameter is required")
    }
    
    // Log pipeline start
    log.info """
    ╔════════════════════════════════════════════════════════════╗
    ║           RNA-Seq Analysis Pipeline v1.0                  ║
    ╠════════════════════════════════════════════════════════════╣
    ║ Input Directory  : ${params.input_dir}
    ║ Reference       : ${params.reference}
    ║ Output Directory: ${params.output_dir}
    ╚════════════════════════════════════════════════════════════╝
    """.stripIndent()
    
    // Run pipeline
    try {
        qc = quality_control(input)
        log.info "Quality control completed successfully"
        
        alignment = align(qc.trimmed)
        log.info "Alignment completed successfully"
    } catch (Exception e) {
        log.error "Pipeline failed: ${e.message}"
        System.exit(1)
    }
    
    // Log completion
    log.info "Pipeline completed successfully!"
    log.info "Results saved to: ${params.output_dir}"
}
```

**Benefits:**
- Clear error messages
- User-friendly logging
- Easy debugging
- Professional appearance

---

## Exercise 7: Complete Production Pipeline (20 minutes)

### Challenge: Build a Professional RNA-Seq Pipeline

**Requirements:**
1. Organize into three modules:
   - `preprocessing.nf` (QC and trimming)
   - `alignment.nf` (alignment and sorting)
   - `counting.nf` (feature counting)

2. Create three subworkflows to wrap them

3. Add parameters:
   - `input_dir`, `reference`, `gtf_file`, `output_dir`
   - `skip_qc`, `align_cpus`, `align_memory`

4. Add error handling and logging

5. Main workflow chains everything

**File structure:**

```
RNA_seq_pipeline/
├── main.nf
├── nextflow.config
├── params.json
└── modules/
    ├── preprocessing.nf
    ├── alignment.nf
    └── counting.nf
```

<details>
<summary>✓ Solution</summary>

**main.nf:**
```groovy
include { preprocessing } from './modules/preprocessing'
include { alignment } from './modules/alignment'
include { counting } from './modules/counting'

params {
    input_dir = "data"
    reference = "reference/hg38.fa"
    gtf_file = "reference/genes.gtf"
    output_dir = "results"
    skip_qc = false
    align_cpus = 8
    align_memory = "16GB"
}

workflow {
    log.info """
    ╔════════════════════════════════════════════════════════════╗
    ║              RNA-Seq Pipeline v2.0                        ║
    ╠════════════════════════════════════════════════════════════╣
    ║ Input Directory  : ${params.input_dir}
    ║ Reference       : ${params.reference}
    ║ GTF File        : ${params.gtf_file}
    ║ Output Dir      : ${params.output_dir}
    ║ Skip QC         : ${params.skip_qc}
    ║ Alignment CPUs  : ${params.align_cpus}
    ╚════════════════════════════════════════════════════════════╝
    """.stripIndent()
    
    // Input validation
    reads = Channel.fromFilePairs("${params.input_dir}/*_{1,2}.fastq.gz")
        .ifEmpty { error("No paired FASTQ files found") }
    
    if (!file(params.reference).exists()) {
        error("Reference not found: ${params.reference}")
    }
    
    if (!file(params.gtf_file).exists()) {
        error("GTF file not found: ${params.gtf_file}")
    }
    
    try {
        // Preprocessing
        if (!params.skip_qc) {
            preprocess_out = preprocessing(reads)
            processed_reads = preprocess_out.trimmed
            log.info "Quality control completed"
        } else {
            processed_reads = reads
            log.info "Skipping quality control"
        }
        
        // Alignment
        reference = Channel.value(params.reference)
        align_out = alignment(processed_reads, reference)
        log.info "Alignment completed"
        
        // Counting
        gtf = Channel.value(params.gtf_file)
        counts = counting(align_out.bam_files, gtf)
        log.info "Feature counting completed"
        
        // Output
        counts.view()
        
    } catch (Exception e) {
        log.error "Pipeline failed: ${e.message}"
        System.exit(1)
    }
    
    log.info "Pipeline completed successfully!"
}
```

**modules/preprocessing.nf:**
```groovy
process FASTQC {
    input: tuple val(id), path(reads)
    output: tuple val(id), path("*.zip")
    script: """
        fastqc ${reads[0]} ${reads[1]}
    """
}

process TRIMMOMATIC {
    input: tuple val(id), path(reads)
    output: tuple val(id), path("*_R{1,2}_paired.fastq.gz")
    script: """
        trimmomatic PE -threads 4 \\
            ${reads[0]} ${reads[1]} \\
            ${id}_R1_paired.fastq.gz ${id}_R1_unpaired.fastq.gz \\
            ${id}_R2_paired.fastq.gz ${id}_R2_unpaired.fastq.gz \\
            ILLUMINACLIP:adapters.fa:2:30:10
    """
}

workflow preprocessing {
    take: reads
    main:
        FASTQC(reads)
        trimmed = TRIMMOMATIC(reads)
    emit:
        trimmed: trimmed
}
```

**modules/alignment.nf:**
```groovy
process ALIGN {
    cpus params.align_cpus
    memory params.align_memory
    
    input:
        tuple val(id), path(reads)
        path reference
    output:
        tuple val(id), path("${id}.bam")
    script: """
        STAR --genomeDir ${reference} \\
             --readFilesIn ${reads[0]} ${reads[1]} \\
             --runThreadN ${task.cpus} \\
             --outSAMtype BAM SortedByCoordinate \\
             --outFileNamePrefix ${id}_
        mv ${id}_Aligned.sortedByCoord.out.bam ${id}.bam
    """
}

process INDEX {
    input:
        tuple val(id), path(bam)
    output:
        tuple val(id), path("${id}.bam"), path("${id}.bam.bai")
    script: """
        samtools index ${bam}
    """
}

workflow alignment {
    take:
        reads
        reference
    main:
        aligned = ALIGN(reads, reference)
        indexed = INDEX(aligned)
    emit:
        bam_files: indexed
}
```

**modules/counting.nf:**
```groovy
process COUNT_FEATURES {
    input:
        tuple val(id), path(bam), path(bai)
        path gtf
    output:
        tuple val(id), path("${id}_counts.txt")
    script: """
        featureCounts -a ${gtf} -o ${id}_counts.txt ${bam}
    """
}

workflow counting {
    take:
        bam_files
        gtf
    main:
        counts = COUNT_FEATURES(bam_files, gtf)
    emit:
        counts: counts
}
```

**nextflow.config:**
```groovy
process {
    publishDir = [path: "${params.output_dir}", mode: 'copy']
    
    withLabel: 'high_mem' {
        memory = '32GB'
        cpus = 16
    }
}

profiles {
    local {
        executor.cpus = 8
        executor.memory = '16GB'
    }
    
    cluster {
        process.executor = 'slurm'
        process.queue = 'normal'
    }
    
    docker {
        docker.enabled = true
        docker.registry = 'quay.io'
    }
}
```

**Run the pipeline:**
```bash
# Basic run
nextflow run main.nf

# Skip QC
nextflow run main.nf --skip_qc true

# Different resources
nextflow run main.nf --align_cpus 16 --align_memory "32GB"

# With profile
nextflow run main.nf -profile cluster

# With parameters file
nextflow run main.nf -params-file params.json
```

</details>

---

## Summary: What You've Practiced

✅ Understanding subworkflows (take/main/emit)  
✅ Creating reusable subworkflows  
✅ Chaining multiple subworkflows  
✅ Using parameters for configuration  
✅ Conditional execution  
✅ Organizing code into modules  
✅ Error handling and validation  
✅ Professional logging  
✅ Building production-ready pipelines  

---

## 🎓 Self-Check

**Can you do these without looking at solutions?**

- [ ] Explain what a subworkflow is
- [ ] Write a subworkflow with take/main/emit
- [ ] Create parameters with defaults
- [ ] Override parameters from command line
- [ ] Use `include` to import modules
- [ ] Add error handling and logging
- [ ] Organize a pipeline into modules
- [ ] Chain multiple subworkflows

**If you can do most of these, you're ready for Week 2!**

---

## 🚀 Congratulations!

You've completed Week 1! You can now:
- ✅ Understand Nextflow concepts (Day 1)
- ✅ Write Groovy code (Day 2)
- ✅ Create processes (Day 3)
- ✅ Build workflows (Day 4)
- ✅ Organize production pipelines (Day 5)

**You're ready to build real bioinformatics pipelines!**

---

## What's Next (Week 2)

- **Day 6:** Practical patterns (splitting, grouping, merging)
- **Day 7:** Real-world example (complete variant calling pipeline)
- **Days 8-14:** Advanced workflows and troubleshooting

---

*Day 5 of 28 - You've completed Week 1 of the Nextflow Mastery Course! 🎉*
