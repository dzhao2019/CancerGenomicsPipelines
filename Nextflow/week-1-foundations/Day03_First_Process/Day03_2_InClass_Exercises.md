# Day 3: Your First Process - Hands-On Exercises

**Time**: 20 minutes for exercises + 15 minutes for solutions  
**Difficulty**: Easy to Moderate  
**Goal**: Write, debug, and understand Nextflow processes

---

## Exercise 1: Understanding Process Structure (5 minutes)

### Part A: Label the Parts

**Given this process, label each part:**

```groovy
process FASTQC_CHECK {
    cpus 2
    memory '4GB'
    
    input:
        path fastq_file
        val sample_id
    
    output:
        path "*.html"
    
    script:
        """
        fastqc ${fastq_file}
        """
}
```

**Identify:**
1. Process name: _______________
2. Directives: _______________
3. Input declaration: _______________
4. Output declaration: _______________
5. Script block: _______________

<details>
<summary>✓ Solution A</summary>

1. **Process name:** `FASTQC_CHECK`
   - Always UPPERCASE by convention
   - Describes what it does

2. **Directives:** `cpus 2` and `memory '4GB'`
   - Tell Nextflow the resource requirements
   - Optional, but good practice

3. **Input declaration:**
   ```groovy
   input:
       path fastq_file
       val sample_id
   ```
   - Declares what the process needs
   - `path` = file, `val` = simple value

4. **Output declaration:**
   ```groovy
   output:
       path "*.html"
   ```
   - Declares what files to capture
   - Pattern matching: all `.html` files

5. **Script block:**
   ```groovy
   script:
       """
       fastqc ${fastq_file}
       """
   ```
   - The actual work
   - Groovy variables interpolated with `${}`

**Key insight:** These 5 parts appear in every process, in this order.

</details>

---

### Part B: Input Types

**Which input type should be used?**

| Input | Type | Why |
|-------|------|-----|
| A FASTQ file path | `path` or `val`? | ___ |
| Sample name "tumor_1" | `path` or `val`? | ___ |
| Number of threads (8) | `path` or `val`? | ___ |
| Reference genome file | `path` or `val`? | ___ |

<details>
<summary>✓ Solution B</summary>

| Input | Type | Why |
|-------|------|-----|
| A FASTQ file path | `path` | It's a file that needs to be processed |
| Sample name "tumor_1" | `val` | It's a string value, not a file |
| Number of threads (8) | `val` | It's a number value, not a file |
| Reference genome file | `path` | It's a file (genome data) |

**Rule of thumb:**
- `path` - Actual files with data (FASTQ, BAM, VCF, etc.)
- `val` - Metadata (strings, numbers, booleans)

</details>

---

## Exercise 2: Fix the Broken Process (5 minutes)

### Part A: Syntax Errors

**Find and fix the errors in this process:**

```groovy
process ALIGN {
    cpus = 8          # Error 1: Wrong syntax
    
    input
        path reads
        val reference   # Error 2: Missing type for reference
    
    output:
        path "*.bam
    
    script:
        """
        bwa mem ${reference} ${reads} > output.bam
        """
}
```

**Identify all errors and rewrite the process correctly.**

<details>
<summary>✓ Solution A</summary>

**Errors found:**

1. **`cpus = 8`** - Should be `cpus 8` (no equals sign)
   - Directives use space, not `=`

2. **`input` missing colon** - Should be `input:`
   - All blocks need colons

3. **`val reference` without `path`** - Missing declaration
   - Each input needs a type

4. **`path "*.bam` - Missing closing quote**
   - Pattern string not properly closed

**Corrected process:**

```groovy
process ALIGN {
    cpus 8
    memory '16GB'
    
    input:
        path reads
        path reference
    
    output:
        path "*.bam"
    
    script:
        """
        bwa mem ${reference} ${reads} > output.bam
        """
}
```

**Key fixes:**
- `cpus 8` (no equals)
- `input:` with colon
- `path reference` (clarified type)
- `path "*.bam"` (closed string)

</details>

---

### Part B: Logic Errors

**The process runs but doesn't work. Why?**

```groovy
process CONVERT_FORMAT {
    input:
        path bam_file
    
    output:
        path "result.vcf"     # Error: Process doesn't create this file!
    
    script:
        """
        samtools view -h ${bam_file} > output.sam
        """
}
```

**What's wrong? How to fix it?**

<details>
<summary>✓ Solution B</summary>

**Problem:**
- Process declares output `path "result.vcf"`
- But script only creates `output.sam`
- Nextflow looks for `result.vcf`, doesn't find it
- Process fails with "No files match pattern"

**Causes:**
1. Wrong file extension in output pattern
2. Script creates different filename than declared output

**Fix (Option 1):** Change output to match what script creates
```groovy
output:
    path "output.sam"  # Match what script creates
```

**Fix (Option 2):** Change script to match output declaration
```groovy
script:
    """
    samtools view -h ${bam_file} > result.vcf
    """
```

**Best practice:** Keep output names meaningful
```groovy
output:
    path "${sample_id}.sam"

script:
    """
    samtools view -h ${bam_file} > ${sample_id}.sam
    """
```

**Key lesson:** Output patterns must match actual files created.

</details>

---

## Exercise 3: Write a Simple Process (10 minutes)

### Part A: Basic Process

**Write a process that:**
- Takes a FASTQ file as input
- Takes a sample ID as input
- Uses FastQC to generate a report
- Outputs the HTML report

**Your process:**

```groovy
process ??? {
    ???
    
    input:
        ???
    
    output:
        ???
    
    script:
        ???
}
```

<details>
<summary>✓ Solution A</summary>

```groovy
process FASTQC {
    cpus 2
    memory '4GB'
    
    input:
        path fastq_file
        val sample_id
    
    output:
        path "*_fastqc.html"
    
    script:
        """
        fastqc ${fastq_file}
        """
}
```

**Explanation:**
- **Name:** `FASTQC` - describes the tool
- **Directives:** Uses 2 CPUs, 4GB RAM (typical for FastQC)
- **Input:** 
  - `path fastq_file` - the sequence file
  - `val sample_id` - metadata (not used in this process, but good practice)
- **Output:** 
  - `path "*_fastqc.html"` - FastQC creates files like `sample_fastqc.html`
- **Script:** Simple call to FastQC

**Alternative (more complete):**
```groovy
process FASTQC {
    cpus 2
    memory '4GB'
    publishDir "results/fastqc"
    
    input:
        path fastq_file
        val sample_id
    
    output:
        path "*_fastqc.html", emit: html
        path "*_fastqc.zip", emit: zip
    
    script:
        """
        fastqc ${fastq_file}
        """
}
```

This version:
- Publishes to a results folder
- Captures both HTML and ZIP outputs
- Gives outputs names (emit) for clarity

</details>

---

### Part B: More Complex Process

**Write a process that:**
- Takes a BAM file and a reference genome as input
- Takes a sample ID as input
- Calls variants using BCFtools
- Names the output VCF with the sample ID
- Outputs the VCF file

**Your process:**

```groovy
process ??? {
    
    input:
        ???
    
    output:
        ???
    
    script:
        ???
}
```

<details>
<summary>✓ Solution B</summary>

```groovy
process CALL_VARIANTS {
    cpus 4
    memory '8GB'
    
    input:
        path bam_file
        path reference
        val sample_id
    
    output:
        path "${sample_id}.vcf"
    
    script:
        """
        bcftools mpileup -f ${reference} ${bam_file} | \
        bcftools call -mv -o ${sample_id}.vcf
        """
}
```

**Explanation:**
- **Inputs:**
  - `path bam_file` - alignment file
  - `path reference` - genome sequence
  - `val sample_id` - for naming output
- **Output:** 
  - `path "${sample_id}.vcf"` - dynamic naming with interpolation
- **Script:**
  - Pipes mpileup → call (memory efficient)
  - `-m` flag for multiallelic calls
  - `-v` flag for variants only
  - `-o` specifies output filename

**Key technique:** Using variable in output pattern:
```groovy
output:
    path "${sample_id}.vcf"  # Groovy interpolation in pattern
```

This creates outputs like:
- `sample1.vcf`
- `sample2.vcf`
- `tumor_001.vcf`

Each one matching the input sample ID!

</details>

---

## Exercise 4: Reading Real Nextflow Code (10 minutes)

### Part A: Understand This Process

**What does this process do?**

```groovy
process ALIGN_AND_SORT {
    cpus 8
    memory '16GB'
    
    input:
        tuple val(sample_id), path(reads), path(reference)
    
    output:
        tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bai")
    
    script:
        """
        bwa mem -t ${task.cpus} ${reference} ${reads} | \
        samtools sort -o ${sample_id}.bam
        
        samtools index ${sample_id}.bam
        """
}
```

**Answer these questions:**
1. What does the process do? (One sentence)
2. What inputs does it have?
3. What outputs does it create?
4. What's the relationship between input and output sample IDs?

<details>
<summary>✓ Solution A</summary>

**1. What does it do?**
Aligns reads to a reference genome, sorts the alignment, and creates an index.

**2. What inputs?**
- `tuple val(sample_id)` - A sample identifier
- `path(reads)` - Sequence reads (FASTQ or similar)
- `path(reference)` - Reference genome

**Note:** Using `tuple` groups these three pieces of data together. They're treated as a unit.

**3. What outputs?**
- `tuple val(sample_id)` - Returns the sample ID (unchanged)
- `path("${sample_id}.bam")` - The aligned, sorted BAM file
- `path("${sample_id}.bai")` - The BAM index file

**4. Relationship between input/output sample IDs?**
The same `sample_id` that comes in is passed through to the output. This keeps the sample ID attached to its data throughout the pipeline.

Example flow:
```
Input: (sample1, reads.fastq, hg38.fa)
    ↓
Process: Align and sort
    ↓
Output: (sample1, sample1.bam, sample1.bai)
```

**Key insight:** By returning the sample ID, the next process knows which sample this BAM came from.

</details>

---

### Part B: What Happens in the Script?

**Trace through the script step-by-step:**

```groovy
script:
    """
    bwa mem -t ${task.cpus} ${reference} ${reads} | \
    samtools sort -o ${sample_id}.bam
    
    samtools index ${sample_id}.bam
    """
```

**With these inputs:**
- `sample_id` = "tumor_001"
- `reads` = "tumor_001.fastq"
- `reference` = "hg38.fa"
- `${task.cpus}` = 8

**What command actually runs?**

<details>
<summary>✓ Solution B</summary>

**Step 1: Variable substitution happens**

Input values are substituted:
- `${task.cpus}` → `8`
- `${reference}` → `hg38.fa`
- `${reads}` → `tumor_001.fastq`
- `${sample_id}` → `tumor_001`

**Step 2: Script becomes:**
```bash
bwa mem -t 8 hg38.fa tumor_001.fastq | \
samtools sort -o tumor_001.bam

samtools index tumor_001.bam
```

**Step 3: What happens:**

Line 1:
```bash
bwa mem -t 8 hg38.fa tumor_001.fastq | samtools sort -o tumor_001.bam
```
- `-t 8` uses 8 CPU cores
- Aligns `tumor_001.fastq` to `hg38.fa`
- Pipes output to `samtools sort`
- Sorted output goes to `tumor_001.bam`

Line 2:
```bash
samtools index tumor_001.bam
```
- Creates index file `tumor_001.bam.bai`

**Files created:**
- `tumor_001.bam` - The alignment
- `tumor_001.bam.bai` - The index

**Key technique:** Using `|` pipe operator is more memory-efficient than creating intermediate files.

</details>

---

## Exercise 5: Debugging Real Errors (10 minutes)

### Part A: Process Fails with "No files match pattern"

**Error message:**
```
ERROR: No files match output pattern: *.vcf
```

**Your process:**
```groovy
process VARIANTS {
    input:
        path bam
    
    output:
        path "*.vcf"
    
    script:
        """
        bcftools call ${bam} > output.txt
        """
}
```

**Debug: Why does it fail? How to fix?**

<details>
<summary>✓ Solution A</summary>

**Why it fails:**
- Process declares it will output `*.vcf` files
- But script only creates `output.txt`
- Nextflow can't find any `.vcf` files
- Process fails

**The fix:**
Either change the output pattern or the script.

**Option 1: Fix the script to create VCF**
```groovy
script:
    """
    bcftools call ${bam} > output.vcf
    """
```

**Option 2: Fix the output pattern**
```groovy
output:
    path "output.txt"
```

**Best practice version:**
```groovy
process VARIANTS {
    input:
        path bam
        val sample_id
    
    output:
        path "${sample_id}.vcf"
    
    script:
        """
        bcftools call ${bam} > ${sample_id}.vcf
        """
}
```

**Lesson:** Always match output pattern to files actually created in script.

</details>

---

### Part B: Process Runs But Output Doesn't Have Sample ID

**Your process works, but you have a problem:**

```groovy
process ALIGN {
    input:
        tuple val(sample_id), path(reads)
    
    output:
        path "output.bam"      # Problem!
    
    script:
        """
        bwa mem reference.fa ${reads} > output.bam
        """
}
```

**Problem:** 
The output is `output.bam` for all samples. Next process won't know which sample this BAM came from!

**Solution:**

```groovy
process ALIGN {
    input:
        tuple val(sample_id), path(reads)
    
    output:
        tuple val(sample_id), path("${sample_id}.bam")
    
    script:
        """
        bwa mem reference.fa ${reads} > ${sample_id}.bam
        """
}
```

**What changed:**
- Output is now a `tuple` with sample ID and BAM file
- BAM file is named with sample ID
- Next process gets sample ID along with the file

**This pattern:**
```groovy
input:
    tuple val(sample_id), ...
output:
    tuple val(sample_id), ...
```

Keeps sample ID traveling through the pipeline!

**Key insight:** Preserve metadata (like sample ID) through processes. It's critical for tracking data.

</details>

---

## Exercise 6: Your First Real Process (15 minutes)

### Challenge: Write a Complete Quality Analysis Process

**Requirements:**
- Name: `QC_REPORT`
- Input: A FASTQ file and sample ID
- Output: An HTML report with sample ID in filename
- Use: FastQC tool
- Make the output publish-ready

**Write the complete process:**

```groovy
process ??? {
    
    ???
    
    input:
        ???
    
    output:
        ???
    
    script:
        ???
}
```

<details>
<summary>✓ Solution</summary>

```groovy
process QC_REPORT {
    cpus 2
    memory '4GB'
    publishDir "results/qc", mode: 'copy'
    
    input:
        path fastq_file
        val sample_id
    
    output:
        path "${sample_id}_fastqc.html"
    
    script:
        """
        fastqc \
            --outdir . \
            --threads ${task.cpus} \
            ${fastq_file}
        
        # FastQC creates file with original name
        # Rename to include sample_id for clarity
        mv *_fastqc.html ${sample_id}_fastqc.html
        """
}
```

**Explanation:**

1. **Process name:** `QC_REPORT` - Clear what it does

2. **Directives:**
   - `cpus 2` - FastQC uses 2 CPUs
   - `memory '4GB'` - Needs 4GB RAM
   - `publishDir "results/qc"` - Copy outputs here
   - `mode: 'copy'` - Copy files (not symlink)

3. **Input:**
   - `path fastq_file` - The sequence file
   - `val sample_id` - For naming output

4. **Output:**
   - `path "${sample_id}_fastqc.html"` - Expected output name

5. **Script:**
   - `--outdir .` - Create output in current directory
   - `--threads ${task.cpus}` - Use assigned CPU cores
   - Rename output to include sample ID

**When called with:**
```
sample_id = "tumor_001"
fastq_file = "tumor_001.fastq"
```

**Creates:**
- `results/qc/tumor_001_fastqc.html` - Published output
- `tumor_001_fastqc.zip` - FastQC zip file (extra)

**Real-world notes:**
- FastQC creates specific filenames
- Renaming ensures consistency
- `publishDir` copies to results folder
- Process is production-ready!

</details>

---

## Summary: What You've Practiced

✅ Understanding process structure (5 parts)  
✅ Fixing syntax errors in processes  
✅ Fixing logic errors (output mismatches)  
✅ Writing simple processes from scratch  
✅ Reading and understanding complex processes  
✅ Debugging real process failures  
✅ Writing production-ready processes  

---

## 🎓 Self-Check

**Can you do these without looking at solutions?**

- [ ] Write a process that takes a file and creates an output
- [ ] Explain the difference between `path` and `val` inputs
- [ ] Fix a process with output pattern mismatch
- [ ] Understand variable substitution in scripts
- [ ] Read a real Nextflow process and explain it
- [ ] Debug a process that's creating the wrong files
- [ ] Use `publishDir` to organize outputs

**If you can do most of these, you're ready for Day 4!**

---

## 🚀 What's Next (Day 4)

Tomorrow you'll learn about **Channels** - How data flows through processes.

You already understand:
- What a process is ✅
- How to write a process ✅
- How processes work ✅

Tomorrow you'll learn:
- How to create channels
- How to connect processes
- How to build complete workflows

Your first workflow is just one day away!

---

*Day 3 of 28 - You've written your first Nextflow process! 🎉*
