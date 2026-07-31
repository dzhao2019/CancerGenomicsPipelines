# Day 3: Nextflow Processes - Quick Reference Card

**Print this for quick lookup while writing processes!**

---

## 🏗️ Process Template (Copy & Modify)

```groovy
process PROCESS_NAME {
    // Directives (optional but recommended)
    cpus 4
    memory '8GB'
    time '2h'
    publishDir "results"
    
    // Input: what you need
    input:
        path input_file
        val sample_id
    
    // Output: what you create
    output:
        path "*.bam"
        path "*.log"
    
    // Script: the work
    script:
        """
        command ${input_file}
        """
}
```

---

## 📋 Process Structure (5 Parts)

### 1. Process Declaration
```groovy
process PROCESS_NAME {
```
- Name is UPPERCASE
- Describes what it does
- Examples: `QUALITY_CHECK`, `ALIGN`, `VARIANTS`

### 2. Directives (Optional)
```groovy
cpus 4              // CPU cores needed
memory '8GB'        // RAM needed
time '2h'          // Max runtime
publishDir "out"   // Copy outputs here
label 'short'      // For grouping
container 'image'  // Docker/Singularity
```

### 3. Input Block
```groovy
input:
    path file          // Input file
    val value          // Simple value
    tuple val(id), path(f)  // Grouped data
    stdin             // From pipe
```

### 4. Output Block
```groovy
output:
    path "*.bam"       // File pattern
    path "report.txt"  // Specific file
    val "success"      // Value
    stdout            // Console output
```

### 5. Script Block
```groovy
script:
    """
    command ${variable}
    """
```

---

## 🔤 Input Types

| Type | What It Is | Example |
|------|-----------|---------|
| `path file` | File to process | `path fastq` |
| `val value` | String/number/boolean | `val sample_id` |
| `tuple` | Grouped data | `tuple val(id), path(f)` |
| `stdin` | From pipe | `stdin` |

---

## 📤 Output Types

| Type | What It Is | Example |
|------|-----------|---------|
| `path "*.bam"` | File(s) matching pattern | All BAM files |
| `path "output.txt"` | Specific file | Exact name |
| `val "value"` | Simple value | String or number |
| `stdout` | Console output | Everything printed |

---

## 📁 Output Glob Patterns

```groovy
"*.bam"             // All .bam files
"*.{bam,vcf}"       // BAM or VCF files
"sample_*.txt"      // Files starting with sample_
"*_aligned.bam"     // Files ending with _aligned.bam
"${sample_id}.bam"  // Dynamic name with interpolation
```

---

## 📊 Common Directives

```groovy
cpus 4                      // How many CPUs
memory '8GB'                // How much RAM
time '2h'                   // Max runtime (h=hours, m=minutes)
publishDir "results"        // Copy outputs to this folder
publishDir "results", mode: 'copy'  // Explicit copy (not link)
label 'short_running'       // Label for grouping
container 'biocontainers/samtools'  // Docker image
```

---

## 🔧 Script Block Basics

```groovy
script:
    """
    # Comments work
    
    # Single command
    fastqc ${fastq_file}
    
    # Multiple commands
    bwa mem ${ref} ${reads} > output.sam
    samtools sort output.sam > output.bam
    
    # Using special variables
    fastqc --threads ${task.cpus} ${fastq}
    
    # Python in script
    python3 << 'EOF'
    print("Hello from Python")
    EOF
    
    # Conditional bash
    if [ -f output.bam ]; then
        echo "BAM file exists"
    fi
    """
```

---

## 🎯 Special Variables

```groovy
${variable}            // Interpolate Groovy variable
${task.cpus}          // Number of CPU cores
${task.memory}        // Memory in GB
${file.baseName}      // Filename without extension
${file.extension}     // File extension
${file.size()}        // File size in bytes
```

---

## 📝 Real Process Examples

### Simple: FastQC
```groovy
process FASTQC {
    cpus 2
    memory '4GB'
    
    input:
        path fastq
    output:
        path "*_fastqc.html"
    script:
        """
        fastqc ${fastq}
        """
}
```

### Intermediate: Alignment
```groovy
process ALIGN {
    cpus 8
    memory '16GB'
    
    input:
        path reads
        path reference
        val sample_id
    output:
        tuple val(sample_id), path("${sample_id}.bam")
    script:
        """
        bwa mem -t ${task.cpus} ${reference} ${reads} | \
        samtools sort -o ${sample_id}.bam
        """
}
```

### Advanced: Multi-step
```groovy
process CALL_VARIANTS {
    cpus 4
    memory '8GB'
    publishDir "results/vcf"
    
    input:
        tuple val(id), path(bam), path(bai)
        path reference
    output:
        tuple val(id), path("${id}.vcf")
    script:
        """
        bcftools mpileup -f ${reference} ${bam} | \
        bcftools call -mv -o ${id}.vcf
        """
}
```

---

## 🐛 Common Errors & Fixes

| Error | Cause | Fix |
|-------|-------|-----|
| `No files match pattern` | Output doesn't match files created | Check script creates right files |
| `script statement is not valid` | Groovy syntax error | Check quotes, brackets |
| `Not enough inputs` | Missing input data | Check workflow passes data |
| `command not found` | Tool not installed | Install tool or use container |
| `permission denied` | Can't read/write files | Check file permissions |

---

## ✅ Process Checklist

Before running your process, check:

- [ ] Process name is UPPERCASE
- [ ] Input block has `:` at end
- [ ] Output block has `:` at end
- [ ] Script block has triple quotes `"""`
- [ ] Output patterns match actual files
- [ ] Groovy variables use `${}`
- [ ] File inputs use `path` type
- [ ] Simple values use `val` type
- [ ] Script ends properly (no missing bracket)

---

## 🔗 How Processes Work with Channels

```groovy
// You provide input
reads = Channel.fromPath("*.fastq")

// Process receives from channel
process ALIGN {
    input:
        path fastq
    output:
        path "*.bam"
    script: """..."""
}

// Run it
aligned = ALIGN(reads)

// Output goes to channel
aligned.view()  // See results
```

---

## 📌 Best Practices

### 1. Name Descriptively
```groovy
process ALIGN_READS        // Good
process PROCESS1           // Bad
```

### 2. Explicit Inputs/Outputs
```groovy
// Good - clear what's needed
input:
    path fastq_file
    path reference
    val sample_id

// Unclear
input:
    path a
    path b
    val c
```

### 3. Dynamic Naming
```groovy
// Good - preserves sample ID
output:
    path "${sample_id}.bam"

// Bad - all samples get same name
output:
    path "output.bam"
```

### 4. Resource Directives
```groovy
// Good - tells Nextflow what's needed
cpus 8
memory '16GB'

// Missing - system guesses wrong
```

### 5. Error Checking
```groovy
// Good - fail fast if tool missing
if ! command -v bwa &> /dev/null; then
    echo "BWA not found!"
    exit 1
fi

// Bad - mysterious error later
bwa mem ...
```

---

## 🎯 Input/Output Patterns

### Pattern: Sequential Data
```groovy
// Input: one file at a time
input:
    path fastq
output:
    path "*.bam"

// Applied to multiple files
// Process runs once per file
```

### Pattern: Paired Data
```groovy
// Input: grouped tuples
input:
    tuple val(id), path(fastq), val(condition)
output:
    tuple val(id), path("${id}.bam")

// Preserves grouping through workflow
```

### Pattern: Multi-file Output
```groovy
input:
    path fastq
output:
    path "*.bam"
    path "*.log"
    path "*.html"

// Creates channel with all matching files
```

---

## 📚 Method Reference for Scripts

### String Methods
```groovy
${file.toString()}      // Full path
${file.baseName}        // Filename without extension
${file.name}            // Filename only
${file.extension}       // File extension
${file.parent}          // Directory
${file.parent.name}     // Directory name
```

### File Properties
```groovy
${file.size()}          // Size in bytes
${file.exists()}        // true/false
${file.isFile()}        // Is it a file?
${file.isDirectory()}   // Is it a directory?
```

---

## 🚀 Quick Start: Three Processes

### Process 1: Input → Processing → Output
```groovy
process QC {
    input: path fastq
    output: path "*.html"
    script: """
        fastqc ${fastq}
    """
}
```

### Process 2: Multiple Inputs → Processing
```groovy
process ALIGN {
    cpus 8
    input:
        path reads
        path reference
    output:
        path "*.bam"
    script: """
        bwa mem -t ${task.cpus} ${reference} ${reads} > output.bam
    """
}
```

### Process 3: Dynamic Output Naming
```groovy
process RENAME {
    input:
        tuple val(id), path(file)
    output:
        path "${id}_final.bam"
    script: """
        mv ${file} ${id}_final.bam
    """
}
```

---

## 💡 Remember These Rules

1. **Process names are UPPERCASE**
2. **All blocks need colons** (input:, output:, script:)
3. **Output patterns must match created files**
4. **Use `${}`** for variable interpolation
5. **Use `path` for files, `val` for values**
6. **Preserve metadata** (like sample_id) through tuples
7. **Directives are optional** (but recommended)

---

## 🎓 By the Numbers

| Concept | Key Point |
|---------|-----------|
| Parts per process | 5 (directives, input, output, script, name) |
| Input types | 4 (path, val, tuple, stdin) |
| Output types | 4 (path, val, stdout, emit names) |
| Glob characters | 3 (`*`, `?`, `[abc]`) |
| Special variables | 6 (task.cpus, task.memory, file.*, etc.) |

---

## 📖 When to Use What

| Situation | Use This |
|-----------|----------|
| One file input | `path file` |
| String/number input | `val value` |
| Multiple related inputs | `tuple val(...), path(...)` |
| Dynamic output name | `path "${id}.bam"` |
| Fixed output name | `path "report.html"` |
| Multiple output files | Multiple `path` lines |
| Need to preserve ID | Return in output tuple |

---

*Print this page and keep it handy while writing processes!*

**Next:** Day 4 - Channels (connecting processes together)
