# Day 3: Your First Nextflow Process

**Learning Time**: 30 minutes  
**Prerequisites**: Days 1-2 completed, understanding of Groovy basics  
**Goal**: Write and understand your first Nextflow processes

---

## 📖 Introduction (3 minutes)

Welcome to Day 3! Today is exciting because you're writing **actual Nextflow code**. By the end of this session, you'll have written multiple working processes that could be used in real bioinformatics workflows.

**What makes today special**: You're shifting from understanding concepts to creating functional pipeline components. Every process you write today is a reusable building block for future workflows.

### What You'll Learn Today

- The anatomy of a Nextflow process
- How to declare inputs and outputs explicitly
- Why processes are different from Python functions
- How to write script blocks with proper interpolation
- Common process patterns for bioinformatics
- How to structure processes for reusability

### What a Process Actually Is

A **process** is a self-contained unit of work that:
- Declares what data it needs (inputs)
- Declares what data it produces (outputs)  
- Defines what computation to perform (script)
- Can run in isolation (different directory, different machine, different container)

Think of it as a specialized function that Nextflow can:
- Execute in parallel across many inputs
- Retry if it fails
- Resume from checkpoints
- Run in containers
- Distribute across clusters

---

## 🎯 Anatomy of a Process (12 minutes)

### 1. The Basic Structure (3 minutes)

Here's the skeleton of every Nextflow process:

```groovy
process processName {
    // Optional: Directives (cpus, memory, container, etc.)
    
    input:
    // Declare what data this process needs
    
    output:
    // Declare what data this process produces
    
    script:
    // The actual command(s) to execute
    """
    shell commands here
    """
}
```

Let's build this up piece by piece with a real example.

### 2. Your Very First Process (3 minutes)

**Goal**: Count the number of reads in a FASTQ file

```groovy
process countReads {
    input:
    path fastq_file
    
    output:
    path "read_count.txt"
    
    script:
    """
    echo "Counting reads in ${fastq_file}"
    grep -c "^@" ${fastq_file} > read_count.txt
    """
}
```

Let's dissect this:

**Process Name**: `countReads`
- Should be descriptive and follow camelCase convention
- Unique within your workflow

**Input Block**:
```groovy
input:
path fastq_file
```
- `path` means this is a file/directory
- `fastq_file` is the variable name you'll use in the script
- When Nextflow runs this process, it will provide the file

**Output Block**:
```groovy
output:
path "read_count.txt"
```
- `path` indicates this is a file output
- `"read_count.txt"` is the filename that will be created
- Nextflow will look for this file after the script runs
- Can use `${variables}` for dynamic names

**Script Block**:
```groovy
script:
"""
echo "Counting reads in ${fastq_file}"
grep -c "^@" ${fastq_file} > read_count.txt
"""
```
- Triple quotes `"""` for multi-line shell script
- `${fastq_file}` interpolates the Groovy variable
- These are normal shell commands
- The file `read_count.txt` must be created (matches output)

### 3. Input Types Explained (2 minutes)

Nextflow has several input qualifiers:

**`val`** - A value (string, number, variable)
```groovy
input:
val sample_id
val num_threads

script:
"""
echo "Processing sample: ${sample_id}"
echo "Using ${num_threads} threads"
"""
```

**`path`** - A file or directory
```groovy
input:
path input_file
path reference_genome
path "input.bam"  // Can specify expected filename

script:
"""
process_tool ${input_file} ${reference_genome}
# Or use the specified name:
samtools view input.bam
"""
```

**`tuple`** - Multiple related items grouped together
```groovy
input:
tuple val(sample_id), path(reads)

script:
"""
# sample_id is the identifier
# reads is the file
fastqc ${reads} --outdir ${sample_id}
"""
```

**`each`** - Combine each item with all inputs (creates combinations)
```groovy
input:
path sample_file
each chromosome

script:
"""
# Runs once for each combination of sample_file × chromosome
process_chromosome ${sample_file} ${chromosome}
"""
```

### 4. Output Types Explained (2 minutes)

Outputs work similarly to inputs:

**`path`** - File(s) produced
```groovy
output:
path "*.bam"                    // Glob pattern - any .bam files
path "${sample_id}.vcf"         // Dynamic filename
path "results/*"                // Directory contents
```

**`val`** - A value to pass forward
```groovy
output:
val sample_id                   // Pass the sample ID forward
val "${sample_id}_processed"    // Computed value
```

**`tuple`** - Multiple outputs grouped
```groovy
output:
tuple val(sample_id), path("*.bam")

// This emits pairs of (sample_id, bam_file)
```

**`stdout`** - Capture standard output
```groovy
output:
stdout

script:
"""
echo "Processing complete"
"""
// The text "Processing complete" becomes the output
```

### 5. A Complete Real-World Example (2 minutes)

**Process**: Run FastQC on a FASTQ file

```groovy
process runFastQC {
    // Directive: specify container (we'll learn more about this later)
    container 'biocontainers/fastqc:v0.11.9_cv8'
    
    // Directive: use 2 CPUs
    cpus 2
    
    input:
    tuple val(sample_id), path(fastq)
    
    output:
    tuple val(sample_id), path("${sample_id}_fastqc.html")
    tuple val(sample_id), path("${sample_id}_fastqc.zip")
    
    script:
    """
    # Create output directory
    mkdir -p fastqc_output
    
    # Run FastQC
    fastqc -t ${task.cpus} -o fastqc_output ${fastq}
    
    # Rename outputs to include sample ID
    mv fastqc_output/*_fastqc.html ${sample_id}_fastqc.html
    mv fastqc_output/*_fastqc.zip ${sample_id}_fastqc.zip
    """
}
```

**New concepts here**:

1. **Directives**: `container` and `cpus` are directives that configure how the process runs
2. **`task.cpus`**: Built-in variable that reflects the cpus directive value
3. **Tuple input/output**: Keeps sample_id paired with its files
4. **Multiple outputs**: Can emit multiple files/values
5. **Shell script complexity**: Real commands with directory management

### 6. Common Patterns You'll Use (2 minutes)

**Pattern 1: Simple transformation**
```groovy
process compressFile {
    input:
    path input_file
    
    output:
    path "${input_file}.gz"
    
    script:
    """
    gzip -c ${input_file} > ${input_file}.gz
    """
}
```

**Pattern 2: With sample tracking**
```groovy
process alignReads {
    input:
    tuple val(sample_id), path(reads)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    script:
    """
    bwa mem ${reference} ${reads} > ${sample_id}.sam
    samtools view -b ${sample_id}.sam > ${sample_id}.bam
    """
}
```

**Pattern 3: Multiple inputs from different sources**
```groovy
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
        -O ${sample_id}.vcf \\
        --dbsnp ${known_sites}
    """
}
```

**Pattern 4: Optional outputs**
```groovy
process analyzeData {
    input:
    path data_file
    
    output:
    path "results.txt"
    path "*.png", optional: true  // May or may not be created
    
    script:
    """
    analyze.sh ${data_file} > results.txt
    # PNG plots are created only for some data types
    """
}
```

---

## 🔗 Python Connection: Functions vs Processes (3 minutes)

Let's compare a Python function with a Nextflow process:

### Python Function Approach

```python
import subprocess
import os

def run_fastqc(fastq_file, output_dir="results"):
    """Run FastQC on a FASTQ file"""
    # Implicit: assumes FastQC is in PATH
    # Implicit: runs in current environment
    # Implicit: uses current directory
    # Implicit: returns... what exactly?
    
    os.makedirs(output_dir, exist_ok=True)
    
    cmd = [
        "fastqc",
        "-o", output_dir,
        fastq_file
    ]
    
    subprocess.run(cmd, check=True)
    
    # Return value is unclear - maybe the output path?
    return os.path.join(output_dir, f"{fastq_file}_fastqc.html")

# Usage
result = run_fastqc("sample1.fastq")
```

**Problems**:
- Runs in the same Python process/environment
- No isolation between runs
- Hard to parallelize safely
- No automatic retry on failure
- Environment dependencies are implicit
- Output tracking is manual
- Resource requirements are unknown

### Nextflow Process Approach

```groovy
process runFastQC {
    container 'biocontainers/fastqc:v0.11.9_cv8'  // Explicit environment
    cpus 2                                         // Explicit resources
    
    input:
    path fastq_file                                // Explicit input
    
    output:
    path "*_fastqc.html"                          // Explicit output
    path "*_fastqc.zip"
    
    script:
    """
    fastqc -t ${task.cpus} ${fastq_file}
    """
}

// Usage in workflow
workflow {
    fastq_ch = Channel.fromPath("*.fastq")
    runFastQC(fastq_ch)  // Automatically parallel for all files
}
```

**Benefits**:
- ✅ Runs in isolated work directory
- ✅ Each execution is independent
- ✅ Automatically parallelized
- ✅ Can retry on failure
- ✅ Environment is explicit (container)
- ✅ Outputs are tracked automatically
- ✅ Resources are declared
- ✅ Can resume failed workflows

### Key Differences Table

| Aspect | Python Function | Nextflow Process |
|--------|----------------|------------------|
| **Execution** | Same process/thread | Isolated directory |
| **Parallelization** | Manual (threading/multiprocessing) | Automatic |
| **Dependencies** | Implicit (system) | Explicit (container) |
| **Inputs** | Function parameters | Declared input types |
| **Outputs** | Return value | Declared output files/values |
| **Failure handling** | Try/except | Automatic retry strategies |
| **Resumability** | Manual checkpointing | Built-in with -resume |
| **Portability** | Depends on environment | Runs anywhere |
| **Resource control** | Manual/implicit | Declared (cpus, memory, time) |

### When to Use Each

**Python Function**:
- Data processing logic
- Statistical computations
- Parsing and transformations
- Algorithm implementation

**Nextflow Process**:
- Running external tools (BWA, GATK, etc.)
- Coordinating multiple tools
- Processing many samples
- Need for parallelization
- Production pipelines

**Best Practice**: Use both!
```groovy
process customAnalysis {
    input:
    path data_file
    
    output:
    path "analysis_results.csv"
    
    script:
    """
    # Call your Python script from the process
    python ${projectDir}/bin/analyze.py ${data_file} > analysis_results.csv
    """
}
```

---

## 💻 Hands-On Exercises (10 minutes)

### Exercise 1: Write Your First Process (3 minutes)

**Task**: Create a process that counts lines in a file and outputs the count.

**Requirements**:
- Process name: `countLines`
- Input: A file (any text file)
- Output: A file called `line_count.txt` containing the count
- Use the `wc -l` command

**Template**:
```groovy
process countLines {
    input:
    // Your input declaration here
    
    output:
    // Your output declaration here
    
    script:
    """
    // Your shell command here
    """
}
```

<details>
<summary>Click to see solution</summary>

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
```

**Explanation**:
- `path input_file`: Declares we expect a file as input
- `path "line_count.txt"`: Declares we'll create this output file
- `wc -l ${input_file}`: Shell command that counts lines
- `> line_count.txt`: Redirects output to the file

**Testing concept** (you'll learn to actually run this later):
If you gave this process a file with 100 lines, `line_count.txt` would contain: `100 filename.txt`
</details>

### Exercise 2: Process with Multiple Inputs (4 minutes)

**Task**: Create a process that aligns sequencing reads to a reference genome.

**Requirements**:
- Process name: `alignReads`
- Input 1: Sample ID (a value, not a file)
- Input 2: FASTQ file (paired with the sample ID)
- Input 3: Reference genome (a file that applies to all samples)
- Output: BAM file named `{sample_id}.bam`
- Use a simplified BWA command

**Template**:
```groovy
process alignReads {
    input:
    // Declare tuple with sample ID and reads
    // Declare reference genome
    
    output:
    // Declare output with dynamic name
    
    script:
    """
    // BWA alignment command
    // Convert SAM to BAM
    """
}
```

**Hints**:
- Use `tuple val(...), path(...)` for the sample and reads together
- Use `path` for the reference
- Create a SAM file first, then convert to BAM
- Use `${sample_id}` in filenames

<details>
<summary>Click to see solution</summary>

```groovy
process alignReads {
    input:
    tuple val(sample_id), path(reads)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    script:
    """
    # Align reads to reference
    bwa mem ${reference} ${reads} > ${sample_id}.sam
    
    # Convert SAM to BAM
    samtools view -b ${sample_id}.sam > ${sample_id}.bam
    
    # Clean up intermediate SAM file
    rm ${sample_id}.sam
    """
}
```

**Explanation**:
- `tuple val(sample_id), path(reads)`: Keeps sample ID with its reads
- `path reference`: Single reference file used for all samples
- Output tuple ensures sample_id stays paired with the BAM
- `${sample_id}.bam`: Dynamic filename based on input
- Shell script chains commands together

**Why this structure?**
- The tuple keeps data associated together
- When this runs on 100 samples, each gets its own unique BAM file
- The reference is the same for all samples (not duplicated)
- Output tuple allows next process to know which BAM belongs to which sample
</details>

### Exercise 3: Process with Multiple Outputs (3 minutes)

**Task**: Create a process that runs FastQC and produces both HTML and ZIP outputs.

**Requirements**:
- Process name: `qualityControl`
- Input: Tuple of sample ID and FASTQ file
- Output 1: HTML report (for viewing)
- Output 2: ZIP file (for further processing)
- Both outputs should include the sample ID in their names

**Template**:
```groovy
process qualityControl {
    input:
    // Tuple with sample ID and FASTQ
    
    output:
    // HTML output with sample ID
    // ZIP output with sample ID
    
    script:
    """
    // Run FastQC
    // Rename outputs to include sample ID
    """
}
```

**Hints**:
- FastQC creates `*_fastqc.html` and `*_fastqc.zip`
- You'll need to rename them to include your sample ID
- Use `mv` to rename files
- Can have multiple `path` outputs

<details>
<summary>Click to see solution</summary>

```groovy
process qualityControl {
    input:
    tuple val(sample_id), path(fastq)
    
    output:
    path "${sample_id}_fastqc.html"
    path "${sample_id}_fastqc.zip"
    
    script:
    """
    # Run FastQC (creates files in current directory)
    fastqc ${fastq}
    
    # FastQC creates files with the input filename as prefix
    # We need to rename them to include our sample ID
    
    # Find the created files and rename them
    mv *_fastqc.html ${sample_id}_fastqc.html
    mv *_fastqc.zip ${sample_id}_fastqc.zip
    """
}
```

**Alternative solution with better control**:
```groovy
process qualityControl {
    input:
    tuple val(sample_id), path(fastq)
    
    output:
    tuple val(sample_id), path("${sample_id}_fastqc.html"), path("${sample_id}_fastqc.zip")
    
    script:
    """
    # Create output directory
    mkdir -p qc_output
    
    # Run FastQC with output directory
    fastqc -o qc_output ${fastq}
    
    # Move and rename outputs
    mv qc_output/*_fastqc.html ${sample_id}_fastqc.html
    mv qc_output/*_fastqc.zip ${sample_id}_fastqc.zip
    """
}
```

**Explanation**:
- Multiple `path` outputs - each on its own line
- Or: Single `tuple` output grouping all outputs
- `${sample_id}` creates unique names for each sample
- Shell script renames FastQC's default output names
- Creating a subdirectory keeps things organized during execution

**Why tuple output?**
If you use tuple output, you keep the sample_id with its results, making it easier to track in downstream processes.
</details>

---

## 🤔 Reflection Activity (4 minutes)

### Question 1: Understanding Process Isolation

Consider this Python function:
```python
processed_samples = []

def process_sample(sample):
    result = analyze(sample)
    processed_samples.append(result)  # Shared state!
    return result
```

**Questions**:
1. What problems could occur if this runs in parallel?
2. How does Nextflow's process isolation solve this?
3. What's the Nextflow way to collect results?

<details>
<summary>Click for answers</summary>

**Answers**:

1. **Parallel problems**:
   - Race conditions on `processed_samples` list
   - Results could be lost or duplicated
   - Need locks/mutexes for thread safety
   - Debugging is difficult

2. **Nextflow isolation**:
   - Each process runs in its own directory
   - No shared state between process executions
   - Results communicated via outputs, not shared variables
   - Parallelization is safe by design

3. **Nextflow collection**:
   ```groovy
   process analyze {
       input:
       val sample
       
       output:
       path "result.txt"
       
       script:
       """
       analyze_tool ${sample} > result.txt
       """
   }
   
   workflow {
       samples = Channel.from(["s1", "s2", "s3"])
       results = analyze(samples)
       all_results = results.collect()  // Collect all outputs
   }
   ```
   
   Results flow through channels, no shared state needed!
</details>

### Question 2: Input/Output Design

You need to create a process that:
- Takes paired-end FASTQ files (read1 and read2)
- Takes a sample ID
- Produces an aligned BAM file

**Design the input and output blocks**:
```groovy
process alignPairedReads {
    input:
    // Your answer here
    
    output:
    // Your answer here
}
```

<details>
<summary>Click for suggested solution</summary>

```groovy
process alignPairedReads {
    input:
    tuple val(sample_id), path(read1), path(read2)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    script:
    """
    bwa mem ${reference} ${read1} ${read2} | \\
        samtools view -b > ${sample_id}.bam
    """
}
```

**Why this design?**
- `tuple` groups related inputs together
- `val(sample_id)` for the identifier
- `path(read1), path(read2)` for the paired files
- `path reference` separate because it's shared across samples
- Output tuple keeps sample_id with the BAM for downstream tracking

**Alternative** - if reads are in a list:
```groovy
input:
tuple val(sample_id), path(reads)  // reads is a list of 2 files
path reference

script:
"""
bwa mem ${reference} ${reads[0]} ${reads[1]} | \\
    samtools view -b > ${sample_id}.bam
"""
```
</details>

### Question 3: Process vs Function

For each scenario, decide: Write as a **Nextflow process** or **Python function**?

**Scenario A**: Calculate GC content percentage from a sequence string
```python
def gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100
```

**Scenario B**: Run BLAST on 1000 protein sequences against a database

**Scenario C**: Parse a VCF file and filter variants by quality score

**Scenario D**: Trim adapters from 500 FASTQ files using Trimmomatic

<details>
<summary>Click for answers</summary>

**Answers**:

**Scenario A: Python function** ✅
- Pure data processing
- No external tools
- Single string input
- Fast computation
- Perfect for Python

**Scenario B: Nextflow process** ✅
```groovy
process runBlast {
    input:
    path protein_seq
    path database
    
    output:
    path "${protein_seq.baseName}.blast"
    
    script:
    """
    blastp -query ${protein_seq} -db ${database} > ${protein_seq.baseName}.blast
    """
}
```
- External tool (BLAST)
- Many inputs (1000 sequences)
- Needs parallelization
- Could take hours serially

**Scenario C: Python function** ✅
```python
def filter_variants(vcf_file, min_quality):
    # Parse and filter logic here
    pass
```
- Data processing/parsing
- Python's strength
- Could be called FROM a Nextflow process though!

**Scenario D: Nextflow process** ✅
```groovy
process trimReads {
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq")
    
    script:
    """
    trimmomatic SE ${reads} ${sample_id}_trimmed.fastq LEADING:3 TRAILING:3
    """
}
```
- External tool (Trimmomatic)
- Many files (500)
- Needs parallelization
- Classic Nextflow use case
</details>

---

## 📝 Key Takeaways

Before moving to Day 4, ensure you understand:

✅ **A process has three main parts**: input, output, script  
✅ **Inputs are explicitly declared** with types (val, path, tuple)  
✅ **Outputs are explicitly declared** and must be created by the script  
✅ **Script blocks use triple quotes** and interpolate with `${}`  
✅ **Processes run in isolation** - each in its own directory  
✅ **Processes are parallelizable** - same process, different inputs, runs simultaneously  
✅ **Tuple keeps related data together** - sample ID with its files  

### The Mental Model: Process as a Factory Station

Think of a process as a **factory workstation**:

- **Input**: Raw materials delivered to the station
- **Script**: The work performed at the station
- **Output**: Finished products leaving the station
- **Isolation**: Each station is independent
- **Parallelization**: Multiple stations of the same type can run simultaneously

```
Input (fastq) → [Process: Quality Check] → Output (html, zip)
Input (fastq) → [Process: Quality Check] → Output (html, zip)  ← Runs in parallel!
Input (fastq) → [Process: Quality Check] → Output (html, zip)  ← Runs in parallel!
```

### Common Process Patterns Summary

**Pattern 1: Simple transformation**
```groovy
process transformFile {
    input: path(file)
    output: path("output.txt")
    script: "tool ${file} > output.txt"
}
```

**Pattern 2: With sample tracking**
```groovy
process processWithID {
    input: tuple val(id), path(file)
    output: tuple val(id), path("${id}_result.txt")
    script: "tool ${file} > ${id}_result.txt"
}
```

**Pattern 3: Multiple inputs**
```groovy
process combineInputs {
    input:
    tuple val(id), path(sample)
    path reference
    output: tuple val(id), path("${id}.out")
    script: "tool -r ${reference} ${sample} > ${id}.out"
}
```

**Pattern 4: Multiple outputs**
```groovy
process multipleOutputs {
    input: path(input)
    output:
    path "result.txt"
    path "report.html"
    path "*.png", optional: true
    script: "complex_tool ${input}"
}
```

---

## 🎯 Ready for Day 4?

Tomorrow, you'll learn about **channels** - the data highways that connect processes. You'll understand how data flows through workflows and enables automatic parallelization.

### Quick Preview

```groovy
// You wrote this today:
process countReads {
    input: path(fastq)
    output: path("count.txt")
    script: "wc -l ${fastq} > count.txt"
}

// Tomorrow you'll learn how to feed it data:
workflow {
    // Channel creates the data stream
    fastq_files = Channel.fromPath("data/*.fastq")
    
    // Process consumes from the channel
    counts = countReads(fastq_files)  // Automatic parallel execution!
}
```

### Process Quick Reference Card

```groovy
// Basic process structure
process myProcess {
    // Optional directives
    cpus 4
    memory '8 GB'
    container 'biocontainers/tool:1.0'
    
    input:
    val simple_value              // A value (string, number)
    path input_file              // A file
    tuple val(id), path(file)    // Related data grouped
    
    output:
    path "output.txt"            // A file created
    val computed_value           // A value
    tuple val(id), path("*.bam") // Multiple things grouped
    stdout                       // Capture printed output
    
    script:
    """
    # Shell commands here
    tool ${input_file} > output.txt
    echo "Processed ${id}"
    """
}
```

---

## ✅ Day 3 Completion Checklist

Before marking Day 3 complete, ensure you can:

- [ ] Write a basic process with input, output, and script
- [ ] Use `path` for file inputs and outputs
- [ ] Use `val` for value inputs
- [ ] Use `tuple` to group related data
- [ ] Interpolate variables in script blocks with `${}`
- [ ] Understand why processes are isolated
- [ ] Explain the difference between a process and a Python function
- [ ] Write processes for real bioinformatics tools

**Self-Test**: Can you write a process that:
- Takes a BAM file and sample ID as input
- Runs `samtools flagstat` on it
- Outputs a stats file named `{sample_id}_stats.txt`

<details>
<summary>Check your answer</summary>

```groovy
process bamStatistics {
    input:
    tuple val(sample_id), path(bam_file)
    
    output:
    tuple val(sample_id), path("${sample_id}_stats.txt")
    
    script:
    """
    samtools flagstat ${bam_file} > ${sample_id}_stats.txt
    """
}
```

If you got this right, you're ready for Day 4! 🎉
</details>

**Completed Day 3?** Update your `PROGRESS.md`! You can now write Nextflow processes! 🚀

**Your progress**: 3/28 days (10.7%) complete

---

## 🔍 Troubleshooting Corner

**Common mistakes when starting**:

1. **Forgetting the colon after input/output/script**
   ```groovy
   input  // ❌ Missing colon
   input: // ✅ Correct
   ```

2. **Not using triple quotes for script**
   ```groovy
   script:
   "single line"      // ✅ OK for one line
   """multi
      line"""         // ✅ Required for multiple lines
   ```

3. **Forgetting ${} in interpolation**
   ```groovy
   script:
   """
   echo $variable     // ⚠️ Shell variable, not Groovy!
   echo ${variable}   // ✅ Groovy variable interpolated
   """
   ```

4. **Output file not created**
   ```groovy
   output:
   path "result.txt"  // Declares this file must exist
   
   script:
   """
   tool input.txt     // ❌ Doesn't create result.txt!
   """
   // This will fail! Always create the declared outputs.
   ```

5. **Using quotes inconsistently**
   ```groovy
   output:
   path "${id}.txt"   // ✅ Correct - interpolation
   path '${id}.txt'   // ❌ No interpolation with single quotes
   ```

---

*Tomorrow: Day 4 - Understanding Channels*

**Excellent progress! See you tomorrow! 🚀**