# Day 3: Your First Nextflow Process

**Learning Time**: 30 minutes  
**Prerequisites**: Completed Days 1-2  
**Goal**: Write, understand, and run your first Nextflow process

---

## 📖 Introduction (2 minutes)

Welcome to Day 3! Over the last two days, you've learned:
- **Day 1:** What Nextflow is and why it matters
- **Day 2:** Groovy syntax to read and write Nextflow code

Today you'll put it all together: **Write your first Nextflow process.**

This is the moment where Nextflow becomes real. You'll go from theory to actually writing working code that processes data.

### What You're About to Do

By the end of this lesson, you'll:
1. Understand process structure completely
2. Write your own process from scratch
3. Run it and see it work
4. Understand what happens at each step

**This is a pivotal moment.** After this, you'll have the power to write real bioinformatics workflows.

---

## 🎯 Learning Objectives

By the end of Day 3, you should be able to:

1. **Explain process anatomy** - Every part and why it exists
2. **Write input declarations** - Tell Nextflow what data you need
3. **Write output declarations** - Tell Nextflow what you're creating
4. **Write script blocks** - Execute real commands
5. **Use variables correctly** - Pass data into your script
6. **Run a process** - Execute it and see results
7. **Debug basic errors** - Understand when/why things go wrong
8. **Understand the flow** - How data moves through a process

---

## 📚 Key Concepts (20 minutes)

### 1. Process Structure: The Anatomy

Every Nextflow process has the same basic structure:

```groovy
process PROCESS_NAME {
    // 1. Directives (optional - instructions to Nextflow)
    cpus 4
    memory '8GB'
    
    // 2. Input (what the process needs)
    input:
        path fastq_file
        val sample_id
    
    // 3. Output (what the process creates)
    output:
        path "*.bam"
    
    // 4. Script (the actual work - can be bash, python, etc.)
    script:
        """
        bwa mem reference.fa ${fastq_file} > ${sample_id}.bam
        """
}
```

**Let's break down each part:**

#### Part 1: Process Declaration
```groovy
process PROCESS_NAME {
```
- `process` is a keyword (always the same)
- `PROCESS_NAME` should be UPPERCASE (convention)
- Name describes what the process does
- Examples: `QUALITY_CHECK`, `ALIGN_READS`, `CALL_VARIANTS`

#### Part 2: Directives (Optional)
```groovy
cpus 4
memory '8GB'
time '2h'
publishDir "results"
label 'short_running'
```

These tell Nextflow how to run the process:
- `cpus` - How many CPU cores it needs
- `memory` - How much RAM it needs
- `time` - Maximum runtime allowed
- `publishDir` - Where to copy output files
- `label` - Tag for grouping similar processes

**You don't always need directives.** Many processes work without them. Nextflow provides defaults.

#### Part 3: Input Block
```groovy
input:
    path fastq_file
    val sample_id
    path reference
```

Declares what data the process needs:
- `path` - Input is a file
- `val` - Input is a simple value (string, number, etc.)
- Names are arbitrary (you choose them)

**Key concept:** These inputs come from channels. When you run the process, Nextflow passes data through channels.

#### Part 4: Output Block
```groovy
output:
    path "*.bam"
    path "*.log"
    val "success"
```

Declares what the process produces:
- `path` - Output file(s) matching a pattern
- `val` - Output value
- Can emit to multiple channels

**Patterns matter:** 
- `"*.bam"` - All `.bam` files
- `"*_aligned.bam"` - Specific naming
- Must match actual files created

#### Part 5: Script Block
```groovy
script:
    """
    bwa mem reference.fa ${fastq_file} > ${sample_id}.bam
    samtools index ${sample_id}.bam
    """
```

The actual work:
- Triple quotes `"""` allow multi-line commands
- Can be bash, Python, R, or any language
- Groovy variables are interpolated with `${}`
- Everything runs in a temporary directory

**Important:** The script is where your actual bioinformatics happens. Nextflow just coordinates it.

---

### 2. Understanding Inputs: Data Types

Nextflow has different input types for different data:

#### Input Type 1: `path` (Files)
```groovy
input:
    path fastq_file    // A file
```
- Receives a file object
- You can use its properties: `${fastq_file.baseName}`, `${fastq_file.size()}`
- Nextflow handles copying the file where it's needed

**Example in script:**
```bash
fastqc ${fastq_file}  # File is copied, then processed
```

#### Input Type 2: `val` (Simple Values)
```groovy
input:
    val sample_id      // A string or number
    val quality_score  // Another value
```
- Receives strings, numbers, booleans
- NOT files (use `path` for files)
- Simple to use: just interpolate

**Example in script:**
```bash
echo "Processing ${sample_id}" > log.txt
threshold=${quality_score}
```

#### Input Type 3: `tuple` (Grouped Data)
```groovy
input:
    tuple val(id), path(fastq), val(condition)
```
- Groups related data together
- `val(id)` - The id is a value
- `path(fastq)` - The fastq is a file
- `val(condition)` - The condition is a value

**Example in script:**
```bash
echo "Sample ${id}: condition=${condition}"
process_fastq.pl ${fastq}
```

#### Input Type 4: `stdin` (Standard Input)
```groovy
input:
    stdin
```
- Receives data on standard input
- Less common in bioinformatics
- Useful for piping between commands

---

### 3. Understanding Outputs: What You Create

Outputs tell Nextflow what to capture and pass to the next process.

#### Output Type 1: `path` (Files)
```groovy
output:
    path "*.bam"           // All .bam files
    path "report.html"     // Specific file
    path "${sample_id}.vcf" // Dynamic name
```

**How it works:**
1. Process creates files in its work directory
2. Nextflow captures files matching the pattern
3. Files become available in a channel
4. Can be passed to next process

**Glob patterns:**
- `*` - Any characters
- `?` - Single character
- `[abc]` - One of a, b, or c
- Examples:
  - `*.bam` - All BAM files
  - `*.{bam,vcf}` - BAM or VCF files
  - `sample_*.fastq` - FASTQ files starting with "sample_"

#### Output Type 2: `val` (Values)
```groovy
output:
    val "success"           // Literal value
    val "${sample_id}"      // Interpolated value
    path "stats.txt"        // But want to capture as value
```

Less common than `path`, but useful for:
- Passing status information
- Count results
- Sample IDs

#### Output Type 3: `stdout` (Standard Output)
```groovy
output:
    stdout
```
- Captures everything printed to console
- Creates a text channel
- Useful for simple tools that print results

**Example:**
```groovy
script:
    """
    echo "Processed: ${sample_id}"
    wc -l ${fastq_file}
    """
```

---

### 4. The Script Block: Where Work Happens

The script block is where bioinformatics happens:

```groovy
script:
    """
    # Comments are allowed (bash comments)
    
    # Step 1: Run a tool
    fastqc ${fastq_file}
    
    # Step 2: Run another tool
    bwa mem ${reference} ${fastq_file} > output.sam
    
    # Step 3: Convert format
    samtools view -b output.sam > output.bam
    
    # Step 4: Index
    samtools index output.bam
    """
```

**Key points:**
- Triple quotes allow multiple lines
- Bash commands (or Python, R, etc.)
- Groovy variables are interpolated with `${}`
- Each command runs in sequence
- If any step fails, process fails

**Multiple languages example:**
```groovy
script:
    """
    # Bash command
    fastqc ${fastq_file}
    
    # Python embedded in bash
    python3 << 'EOF'
    import sys
    sample = "${sample_id}"
    print(f"Processing {sample}")
    EOF
    """
```

---

### 5. Example: Simple Quality Control Process

Let's build a real example step by step:

```groovy
// Step 1: Process name describes what it does
process QUALITY_CHECK {
    
    // Step 2: Directives (optional, but good practice)
    cpus 2
    memory '4GB'
    
    // Step 3: Input - what we need
    input:
        path fastq_file
        val sample_id
    
    // Step 4: Output - what we create
    output:
        path "${sample_id}_fastqc.html"
    
    // Step 5: The script - do the work
    script:
        """
        fastqc \
            --outdir . \
            --threads ${task.cpus} \
            ${fastq_file}
        
        # Rename output to match our sample_id
        mv *_fastqc.html ${sample_id}_fastqc.html
        """
}
```

**Breaking it down:**

1. **Name:** `QUALITY_CHECK` - Clear what it does
2. **Directives:** Uses 2 CPUs, needs 4GB RAM
3. **Input:** Takes a FASTQ file and a sample ID
4. **Output:** Produces an HTML report with the sample ID in the filename
5. **Script:** 
   - Runs FastQC with proper options
   - Uses `${task.cpus}` (Nextflow variable for CPU count)
   - Renames output to match our naming scheme

---

### 6. Process Communication: How Data Flows

Processes don't directly talk to each other. They communicate through channels:

```groovy
// Workflow (you'll learn this on Day 4, but shown for context)
workflow {
    // Create input channel with FASTQ files
    reads = Channel.fromPath("data/*.fastq")
    
    // Run process for each file
    qc_results = QUALITY_CHECK(reads, reads.map { it.baseName })
    
    // Results are now in qc_results channel
    // Can be passed to next process
}
```

**The flow:**
1. Input channel has multiple FASTQ files
2. For each file, run QUALITY_CHECK
3. QUALITY_CHECK creates output files
4. Output files go to output channel
5. Next process receives from that channel

**Key insight:** Processes don't wait for all inputs. They start processing as soon as data arrives.

---

### 7. Real Example: Alignment Process

A more complex example with multiple tools:

```groovy
process ALIGN_READS {
    cpus 8
    memory '16GB'
    time '4h'
    publishDir "results/bam"
    
    input:
        path reads
        path reference
        val sample_id
    
    output:
        tuple val(sample_id), path("*.bam"), path("*.bai")
    
    script:
        """
        # Align reads to reference
        bwa mem \
            -t ${task.cpus} \
            -M \
            ${reference} \
            ${reads} \
            | samtools view -b \
            > ${sample_id}.bam
        
        # Sort BAM file
        samtools sort \
            -o ${sample_id}.sorted.bam \
            ${sample_id}.bam
        
        # Index the BAM file
        samtools index ${sample_id}.sorted.bam
        
        # Clean up unsorted BAM
        rm ${sample_id}.bam
        
        # Rename sorted to final
        mv ${sample_id}.sorted.bam ${sample_id}.bam
        mv ${sample_id}.sorted.bam.bai ${sample_id}.bam.bai
        """
}
```

**What's happening:**
1. Takes reads, reference, and sample ID
2. Uses multiple CPU cores for alignment
3. Pipes output through samtools (memory efficient)
4. Sorts and indexes the BAM file
5. Cleans up intermediate files
6. Outputs the final BAM and index
7. `publishDir` copies results to a results folder

---

### 8. Common Process Patterns

#### Pattern 1: Quality Control
```groovy
process QC {
    input:
        path reads
    output:
        path "*.html"
    script:
        """
        fastqc ${reads}
        """
}
```

#### Pattern 2: Transformation
```groovy
process FORMAT_CONVERSION {
    input:
        path sam_file
    output:
        path "*.bam"
    script:
        """
        samtools view -b ${sam_file} > output.bam
        """
}
```

#### Pattern 3: Analysis
```groovy
process STATISTICS {
    input:
        path bam_file
    output:
        path "stats.txt"
    script:
        """
        samtools stats ${bam_file} > stats.txt
        """
}
```

#### Pattern 4: Aggregation (Multiple Inputs)
```groovy
process MERGE_RESULTS {
    input:
        path bam_files
    output:
        path "merged.bam"
    script:
        """
        samtools merge merged.bam ${bam_files}
        """
}
```

---

### 9. Understanding the Work Directory

When a process runs, Nextflow creates a temporary directory:

```
work/
├── ab/
│   └── 1234567890abcdef1234567890abcdef/
│       ├── input_files/          (links to input files)
│       ├── script.sh             (the actual script)
│       ├── output_file.bam       (generated output)
│       └── .command.log          (stdout/stderr)
```

**Important facts:**
- Each process run gets its own work directory
- Inputs are linked (not copied) for efficiency
- Outputs stay there until Nextflow collects them
- Logs are saved for debugging
- `-resume` reuses work directories

**This is why:**
- Processes are portable (everything is self-contained)
- They're parallelizable (separate directories, no conflicts)
- Resumable (if you run again, Nextflow knows what succeeded)

---

### 10. Debugging a Process: Common Errors

#### Error 1: File Not Found
```
ERROR: No files match pattern: *.bam
```
**Cause:** Process created output but it doesn't match the pattern  
**Fix:** Check script - does it create that file?

#### Error 2: Script Syntax Error
```
ERROR: script statement is not valid
```
**Cause:** Groovy syntax error in script block  
**Fix:** Check for missing quotes, brackets

#### Error 3: Input Not Provided
```
ERROR: Not enough inputs provided
```
**Cause:** Process has input declarations but workflow doesn't provide them  
**Fix:** Check workflow - is it passing correct data?

#### Error 4: File Permissions
```
ERROR: permission denied
```
**Cause:** Process can't read/write files  
**Fix:** Check file permissions, directories

#### Error 5: Tool Not Found
```
command not found: bwa
```
**Cause:** Tool isn't installed  
**Fix:** Install tool, use container (Day 13)

---

## 🔗 Connecting to Python

### How a Process is Like a Python Function

```python
# Python function
def quality_check(fastq_file, sample_id):
    """Run quality check on a FASTQ file."""
    command = f"fastqc {fastq_file}"
    os.system(command)
    return "report.html"

# Call it
result = quality_check("sample.fastq", "s1")
```

### Nextflow Process Equivalent
```groovy
process QUALITY_CHECK {
    input:
        path fastq_file
        val sample_id
    output:
        path "report.html"
    script:
        """
        fastqc ${fastq_file}
        """
}
```

**Similarities:**
- Input parameters
- Processing step
- Output result
- Can be called with data

**Differences:**
- Process is declarative (describes structure)
- Input/output are explicit
- Runs external tools (not Python code)
- Automatically parallelizes

---

## ✅ Completion Checklist

- [ ] I understand the 5 parts of a process
- [ ] I can explain what each directive does
- [ ] I know the difference between `path` and `val` inputs
- [ ] I understand output patterns (`*.bam`)
- [ ] I can write a simple script block
- [ ] I know what the work directory is
- [ ] I can identify common process errors
- [ ] I feel ready to write my first process
- [ ] I understand how processes communicate

---

## 🔑 Key Takeaways

### The Process Structure (5 Parts)
1. **Directives** - Tell Nextflow how to run it
2. **Input** - Declare what you need
3. **Output** - Declare what you create
4. **Script** - Do the actual work
5. **Process Name** - Describe what it does

### Most Important Rules
1. Process names are UPPERCASE
2. Inputs tell Nextflow what data you need
3. Outputs tell Nextflow what files matter
4. Script is where the work happens
5. Use `${}` to interpolate variables

### Remember
- Processes are isolated units
- They're reusable and parallelizable
- They can run anywhere (laptop, cluster, cloud)
- The same process works for 1 file or 1 million

---

## 🚀 Ready for Exercises?

You now understand process structure completely. Time to write your first one!

The exercises will have you:
1. Fix broken processes
2. Write simple processes
3. Understand real processes
4. Debug actual errors

Let's do it!

---

*This is Day 3 of 28. You're building your first real Nextflow code!*
