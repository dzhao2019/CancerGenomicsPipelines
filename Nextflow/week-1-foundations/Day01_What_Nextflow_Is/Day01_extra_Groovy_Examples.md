# Groovy Code Examples for Day 1

These examples show Groovy code that you'll see in Nextflow. You don't need to understand every detail yet—this is just to familiarize your eyes with the syntax before Day 2's deep dive.

---

## Example 1: Simple Process (The Foundation)

This is the simplest possible Nextflow process. Even if you don't understand every line, you'll recognize the structure in every Nextflow workflow you write:

```groovy
// This is a Nextflow process written in Groovy
process HELLO_WORLD {
    // Output: where does the result go?
    output:
        stdout
    
    // Script: what actually runs
    script:
        """
        echo "Hello, Nextflow!"
        """
}

// The workflow: how to use the process
workflow {
    HELLO_WORLD()
}
```

**Key observations:**
- `process` defines a computational task
- `output: stdout` means send console output as the result
- `script:` contains shell commands (bash, not Groovy!)
- Triple quotes `"""` allow multi-line strings
- `workflow {}` orchestrates processes

---

## Example 2: Process with Input (Processing Data)

This is more realistic. The process takes a filename as input:

```groovy
// Process that takes a FASTQ file as input
process COUNT_LINES {
    // Input: what does this process need?
    input:
        path input_file
    
    // Output: what does it produce?
    output:
        path "line_count.txt"
    
    // Script: the work
    script:
        """
        wc -l ${input_file} > line_count.txt
        """
}

// Use the process
workflow {
    // Create input data
    my_file = file("data/sample.fastq")
    
    // Run the process
    COUNT_LINES(my_file)
}
```

**New concepts:**
- `path input_file`: Input is a file
- `${input_file}`: Groovy string interpolation (like Python f-strings)
- `path "line_count.txt"`: Output is a file
- `file()`: Creates a file object

---

## Example 3: Processing Multiple Files (The Power of Channels)

This is where Nextflow gets powerful. The same code processes 1 file or 1000 files automatically:

```groovy
process PROCESS_FASTQ {
    input:
        path fastq_file
    
    output:
        path "*.stats"
    
    script:
        """
        # Get filename without extension
        SAMPLE=\$(basename ${fastq_file} .fastq)
        
        # Count lines and sequences
        NUM_READS=\$(($(wc -l < ${fastq_file}) / 4))
        
        # Save stats
        echo "Sample: \$SAMPLE" > \$SAMPLE.stats
        echo "Reads: \$NUM_READS" >> \$SAMPLE.stats
        """
}

workflow {
    // Create a channel with multiple files
    // This is the magic—Nextflow will process each file in parallel!
    fastq_files = Channel.fromPath("data/*.fastq")
    
    // Run the process
    // If 10 files: runs 10 times in parallel
    // If 1000 files: schedules them smartly
    PROCESS_FASTQ(fastq_files)
}
```

**Key insight:**
- `Channel.fromPath()` creates a stream of files
- Passing this to a process makes it run multiple times automatically
- This is 100x faster than Python's `for` loop

---

## Example 4: Connecting Two Processes (Building Pipelines)

Watch how one process's output becomes the next process's input:

```groovy
// Process 1: Quality control
process FASTQC {
    input:
        path fastq
    
    output:
        path "*.html"
    
    script:
        """
        fastqc ${fastq}
        """
}

// Process 2: Generate a summary
process SUMMARIZE {
    input:
        path html_file
    
    output:
        path "summary.txt"
    
    script:
        """
        echo "Processed: ${html_file}" > summary.txt
        """
}

// Workflow: Chain them together
workflow {
    // Input files
    samples = Channel.fromPath("data/*.fastq")
    
    // Process 1: Run QC on all samples
    qc_results = FASTQC(samples)
    
    // Process 2: Summarize results
    // This automatically waits for FASTQC to complete!
    SUMMARIZE(qc_results)
}
```

**Key concept:**
- `qc_results = FASTQC(samples)` captures the output as a channel
- `SUMMARIZE(qc_results)` uses that output
- Nextflow automatically manages dependencies

---

## Example 5: Using Tuples (Keeping Related Data Together)

Sometimes you need to process pairs of files (like paired-end reads). Groovy tuples do this:

```groovy
process ALIGN_PAIRED_READS {
    input:
        // Tuple: keeps sample_id and read files together
        tuple val(sample_id), path(read1), path(read2)
    
    output:
        // Output is also a tuple
        tuple val(sample_id), path("*.bam")
    
    script:
        """
        bwa mem reference.fa ${read1} ${read2} | samtools view -b > ${sample_id}.bam
        """
}

workflow {
    // Create tuples: [sample_id, read1, read2]
    sample_pairs = Channel.fromFilePairs("data/*_{1,2}.fastq.gz")
    
    // Process each pair
    // Sample ID stays with the BAM file
    ALIGN_PAIRED_READS(sample_pairs)
}
```

**Why tuples matter:**
- Keep related data together
- Sample ID travels with the data through your pipeline
- Prevents mix-ups when processing many files

---

## Example 6: Parameters (Making Workflows Configurable)

Real workflows take parameters. This allows one workflow to work in many scenarios:

```groovy
// Define parameters with defaults
params {
    input_dir = "data/"
    output_dir = "results/"
    quality_threshold = 30
}

process FILTER_READS {
    input:
        path fastq
    
    output:
        path "*.filtered.fastq"
    
    script:
        // Use the parameter value
        """
        echo "Filtering with quality threshold: ${params.quality_threshold}"
        
        # Some filtering tool that uses the parameter
        seqtk seq -q ${params.quality_threshold} ${fastq} > output.filtered.fastq
        """
}

workflow {
    // Use parameter to find input files
    samples = Channel.fromPath("${params.input_dir}/*.fastq")
    
    // Run the process
    filtered = FILTER_READS(samples)
    
    // Optionally publish results
    filtered.view()
}
```

**This is powerful because:**
- Change behavior without editing code: `nextflow run pipeline.nf --quality_threshold 20`
- Same workflow works for different data
- Easy to share with colleagues

---

## Example 7: Using Groovy Features (Lists and Maps)

Groovy has Python-like data structures:

```groovy
// Lists (like Python lists)
tools = ["fastqc", "bwa", "samtools", "bcftools"]

// Groovy can iterate like Python
tools.each { tool ->
    println "Tool: ${tool}"
}

// Maps (like Python dictionaries)
config = [
    genome: "hg38",
    threads: 8,
    memory: "16GB"
]

// Access like Python dict
println "Genome: ${config.genome}"
println "Threads: ${config.threads}"

// Can also use dot notation
println "Memory: ${config.memory}"

// Combining: List of maps (common in Nextflow)
samples = [
    [id: "sample1", fastq: "sample1.fastq", condition: "control"],
    [id: "sample2", fastq: "sample2.fastq", condition: "treatment"],
]

samples.each { sample ->
    println "Processing ${sample.id} (${sample.condition})"
}
```

**Groovy feels familiar:**
- Similar to Python lists and dicts
- Slightly different syntax, same concepts
- Very natural for data handling

---

## Example 8: String Interpolation (Use This CONSTANTLY)

String interpolation with `${}` is one of the most important Groovy features you'll use:

```groovy
// Basic interpolation
sample_id = "sample_123"
fastq_file = "data/${sample_id}.fastq"
println "Processing: ${fastq_file}"

// Property access in interpolation
input_file = file("data/sample.fastq")
println "File: ${input_file}"
println "File name: ${input_file.baseName}"  // without extension
println "File size: ${input_file.size()}"     // in bytes

// Expressions in interpolation
reads = 1000000
println "Processed ${reads / 1000000} million reads"

// Very common in Nextflow script blocks
process COUNT_READS {
    input:
        path fastq
    
    script:
        """
        # Groovy interpolation happens BEFORE the script runs
        # This is already substituted with actual values
        echo "Counting reads in ${fastq}"
        
        NUM_READS=\$(($(wc -l < ${fastq}) / 4))
        echo "Total reads: \$NUM_READS"
        """
}
```

**Critical concept:**
- Groovy `${}` interpolation happens in Nextflow
- Shell `${}` (escaped as `\${}`) happens in bash
- Watch the escaping carefully!

---

## Example 9: Closures (Advanced, But Important)

Closures are Groovy's lambda functions. You'll see them in Nextflow channel operations:

```groovy
// Simple closure: a block of code
my_list = [1, 2, 3, 4, 5]

// Transform each element
doubled = my_list.collect { it * 2 }
println doubled  // [2, 4, 6, 8, 10]

// Filter elements
evens = my_list.findAll { it % 2 == 0 }
println evens   // [2, 4]

// In Nextflow, you'll use closures with channels
fastq_files = Channel.fromPath("data/*.fastq")

// Transform: remove extension
names = fastq_files.map { file ->
    file.baseName
}

// Filter: only large files
large_files = fastq_files.filter { file ->
    file.size() > 1000000
}

// Chain operations
processed = fastq_files
    .map { file -> [file.baseName, file] }      // Create [name, file] tuple
    .filter { name, file -> file.size() > 1000 } // Keep large files
    .map { name, file -> name }                  // Keep just the name
```

**Closures will make more sense later, but:**
- `{ it * 2 }` is a closure that doubles its input
- `{ file -> file.size() }` is a closure that accesses the file size
- Very powerful for data transformations

---

## Example 10: A Complete Realistic Pipeline

Here's a mini pipeline that shows many concepts together:

```groovy
// Parameters
params {
    input_dir = "raw_data/"
    output_dir = "results/"
    sample_pattern = "*.fastq.gz"
}

// Process 1: Quality control
process QC {
    publishDir "${params.output_dir}/qc"
    
    input:
        tuple val(sample_id), path(fastq)
    
    output:
        tuple val(sample_id), path("*.html")
    
    script:
        """
        fastqc ${fastq}
        """
}

// Process 2: Alignment
process ALIGN {
    input:
        tuple val(sample_id), path(fastq)
        path genome
    
    output:
        tuple val(sample_id), path("*.bam")
    
    script:
        """
        bwa mem ${genome} ${fastq} | samtools view -b > ${sample_id}.bam
        """
}

// Process 3: Merge results
process MERGE_RESULTS {
    publishDir "${params.output_dir}"
    
    input:
        path bam_files
    
    output:
        path "all_samples.txt"
    
    script:
        """
        echo "Processed files:" > all_samples.txt
        ls -lh *.bam >> all_samples.txt
        """
}

// Workflow
workflow {
    // Input: FASTQ files
    samples = Channel.fromFilePairs("${params.input_dir}/${params.sample_pattern}")
    
    // Reference genome
    genome = file("reference/hg38.fa")
    
    // Step 1: QC (parallelized across all samples)
    qc_results = QC(samples)
    
    // Step 2: Alignment (parallelized across all samples)
    aligned = ALIGN(samples, genome)
    
    // Step 3: Merge (waits for alignment to complete)
    MERGE_RESULTS(aligned.map { id, bam -> bam }.collect())
}
```

**This shows:**
- Parameters (default values)
- Multiple processes
- Input/output handling
- Data flow between processes
- Automatic parallelization
- Publishing results

---

## Groovy Syntax Cheat Sheet (Quick Reference)

```groovy
// Variables (no type declaration needed)
name = "sample"
count = 100
is_valid = true

// Strings
single_quoted = 'literal string'
double_quoted = "can use ${variables}"
multi_line = """
This spans
multiple
lines
"""

// Collections
list = [1, 2, 3, 4, 5]
map = [name: "Alice", age: 30]

// Accessing
first_item = list[0]
name_value = map["name"]
name_value = map.name  // same as above

// Iteration
list.each { item ->
    println item
}

// Transformation
doubled = list.map { it * 2 }
evens = list.findAll { it % 2 == 0 }

// File objects
file = file("data/sample.fastq")
file.baseName          // filename without extension
file.size()            // file size in bytes
file.exists()          // boolean

// Interpolation
message = "Hello, ${name}!"
path = "results/${sample_id}/output.txt"

// Control flow (similar to Python)
if (count > 10) {
    println "Large count"
} else {
    println "Small count"
}

for (i = 0; i < 5; i++) {
    println i
}
```

---

## Next Steps

You've now seen Groovy code examples. Tomorrow (Day 2), you'll:
- Learn Groovy syntax in depth
- Understand string interpolation completely
- Learn about closures and their power
- Compare everything to Python

But for now, just recognize that:
- Groovy is similar to Python but with different syntax
- The patterns will become familiar quickly
- You don't need to be a Groovy expert
- These syntax patterns appear in every Nextflow workflow

Keep these examples nearby for reference!
