# Day 1 Hands-On Exercises: Solution Guide

This document contains the detailed solutions and explanations for all Day 1 exercises.

---

## Exercise 1: Identify the Right Tool

### Complete Solutions

#### Scenario A: Data Transformation

**Your Answer Should Be:**
```
Tool: PYTHON
Reason: This is pure data processing and statistical work
```

**Detailed Explanation:**

This task involves:
1. **Data filtering** - Select genes based on expression threshold
2. **Normalization** - Mathematical transformation (CPM)
3. **Visualization** - Creating plots
4. **Single input, single output** - One file in, one file out

**Why Python is better:**
- Libraries like `pandas` excel at data manipulation
- `matplotlib` makes beautiful visualizations
- Statistical operations are straightforward
- Single file processing (no parallelization needed)

**Example Implementation:**
```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the expression matrix
expr_data = pd.read_csv('gene_expression.csv', index_col=0)

# Step 1: Filter low-expression genes (< 5 counts in any sample)
# Keep genes with max expression >= 5
expr_filtered = expr_data[expr_data.max(axis=1) >= 5]

# Step 2: Normalize (counts per million)
expr_normalized = expr_filtered.div(expr_filtered.sum(axis=0)) * 1e6

# Step 3: Create visualizations
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Expression distribution
axes[0, 0].hist(expr_normalized.sum(axis=1), bins=50)
axes[0, 0].set_xlabel('Gene Total Expression (CPM)')
axes[0, 0].set_ylabel('Frequency')
axes[0, 0].set_title('Distribution of Gene Expression')

# Plot 2: Sample coverage
axes[0, 1].bar(expr_normalized.columns, expr_normalized.sum(axis=0))
axes[0, 1].set_xlabel('Sample')
axes[0, 1].set_ylabel('Total Counts')
axes[0, 1].set_title('Coverage Per Sample')

# Plot 3: Heatmap of top genes
top_genes = expr_normalized.nlargest(20, 'sample_1')  # Example
im = axes[1, 0].imshow(top_genes.values, aspect='auto', cmap='RdYlBu_r')
axes[1, 0].set_title('Top 20 Genes')
axes[1, 0].set_xlabel('Samples')
axes[1, 0].set_ylabel('Genes')
plt.colorbar(im, ax=axes[1, 0])

# Plot 4: PCA-like scatter (simplified)
axes[1, 1].scatter(expr_normalized.iloc[:, 0], expr_normalized.iloc[:, 1])
axes[1, 1].set_xlabel('Sample 1 Expression')
axes[1, 1].set_ylabel('Sample 2 Expression')
axes[1, 1].set_title('Sample Correlation')

plt.tight_layout()
plt.savefig('expression_analysis.png', dpi=300)
print(f"Filtered to {len(expr_filtered)} genes from {len(expr_data)}")
print("Visualizations saved!")

# Step 4: Save cleaned data
expr_normalized.to_csv('cleaned_expression.csv')
```

**When Python is best:**
- ✅ Data science/statistics
- ✅ Complex transformations
- ✅ Visualization
- ✅ Single file or in-memory processing
- ✅ Need for mathematical operations

**If it were Nextflow:**
If you had 500 samples to process with the same analysis, then Nextflow would help you parallelize. But for one file, Python is simpler and more appropriate.

---

#### Scenario B: Variant Calling Pipeline

**Your Answer Should Be:**
```
Tool: NEXTFLOW
Reason: Multiple external tools, many samples, need parallelization
```

**Detailed Explanation:**

This task involves:
1. **External tools** - FastQC, BWA, bcftools (not Python-based)
2. **Multiple samples** - 500 FASTQ files
3. **Sequential dependencies** - Each step depends on the previous
4. **Parallelization opportunity** - Can run all samples' step 1 in parallel before step 2

**Why Nextflow is better:**
- Designed to coordinate external tools
- Automatic parallelization across 500 samples
- Built-in resumability (if step 3 fails at sample 400, resume from there)
- Scales to cluster/cloud with config change only
- Each process can run with different software versions

**Example Nextflow Workflow:**

```groovy
// Define parameters
params {
    input_dir = "fastq_files/"
    output_dir = "variants/"
    genome = "reference/hg38.fa"
    threads = 8
}

// Process 1: Quality control
process FASTQC {
    input:
        tuple val(sample_id), path(fastq)
    output:
        tuple val(sample_id), path("*.html")
    script:
        """
        fastqc --threads ${params.threads} ${fastq}
        """
}

// Process 2: Alignment
process ALIGN {
    input:
        tuple val(sample_id), path(fastq)
        path reference
    output:
        tuple val(sample_id), path("*.bam")
    script:
        """
        bwa mem -t ${params.threads} ${reference} ${fastq} | \
        samtools view -b -h -o ${sample_id}.bam
        """
}

// Process 3: Sort BAM
process SORT_BAM {
    input:
        tuple val(sample_id), path(bam)
    output:
        tuple val(sample_id), path("*.sorted.bam")
    script:
        """
        samtools sort -o ${sample_id}.sorted.bam ${bam}
        samtools index ${sample_id}.sorted.bam
        """
}

// Process 4: Variant calling
process CALL_VARIANTS {
    input:
        tuple val(sample_id), path(bam)
        path reference
    output:
        tuple val(sample_id), path("*.vcf")
    script:
        """
        bcftools mpileup -f ${reference} ${bam} | \
        bcftools call -m -o ${sample_id}.vcf
        """
}

// Process 5: Filter variants
process FILTER_VARIANTS {
    publishDir "${params.output_dir}/${sample_id}", mode: 'copy'
    
    input:
        tuple val(sample_id), path(vcf)
    output:
        tuple val(sample_id), path("*.filtered.vcf")
    script:
        """
        bcftools filter -i 'QUAL>30 && DP>10' \
        ${vcf} -o ${sample_id}.filtered.vcf
        """
}

// Process 6: Summary report
process SUMMARY {
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
        path all_vcfs
    output:
        path "summary.txt"
    script:
        """
        echo "Analysis Summary" > summary.txt
        echo "===============" >> summary.txt
        echo "Total samples processed: \$(echo '${all_vcfs}' | wc -w)" >> summary.txt
        echo "Total variants called:" >> summary.txt
        bcftools merge ${all_vcfs} | grep -v '^#' | wc -l >> summary.txt
        """
}

// Workflow
workflow {
    // Load FASTQ files
    samples = Channel.fromFilePairs("${params.input_dir}/*_{1,2}.fastq.gz")
    
    // Load reference
    reference = file(params.genome)
    
    // Execute pipeline
    qc = FASTQC(samples)
    aligned = ALIGN(samples, reference)
    sorted = SORT_BAM(aligned)
    variants = CALL_VARIANTS(sorted, reference)
    filtered = FILTER_VARIANTS(variants)
    
    // Final summary (waits for all samples)
    SUMMARY(filtered.map { id, vcf -> vcf }.collect())
}
```

**Execution behavior:**
```
Sample 1: QC → ALIGN → SORT → VARIANTS → FILTER
Sample 2: QC → ALIGN → SORT → VARIANTS → FILTER  (parallel to sample 1)
Sample 3: QC → ALIGN → SORT → VARIANTS → FILTER  (parallel to sample 1-2)
...
Sample 500: QC → ALIGN → SORT → VARIANTS → FILTER

SUMMARY (waits for all to complete)
```

**Time comparison:**
- Sequential (Python): 500 × (5 + 30 + 5 + 15 + 10) = 32,500 minutes = 22 days
- Parallel (Nextflow): ~55 minutes (assuming 10 parallel samples)
- **Speedup: ~600x faster!**

**When Nextflow is best:**
- ✅ Multiple external tools (FastQC, BWA, SAMtools, bcftools)
- ✅ Many samples (parallelization benefit)
- ✅ Need resumability
- ✅ Multiple analysis steps
- ✅ Output organization
- ✅ Want to scale to cluster/cloud

---

#### Scenario C: Processing One Large File

**Your Answer Should Be:**
```
Tool: PYTHON (or combination)
Reason: Single file, data science operations, complex logic
```

**Detailed Explanation:**

This task involves:
1. **Single file input** - One 200M read BAM file
2. **Statistical analysis** - Calculating complex metrics
3. **Visualization needed** - Coverage plots, junction diagrams
4. **Data science operations** - Not just tool coordination
5. **Complex logic** - Finding junctions, aggregating by chromosome

**Why Python is better:**
- One sample means no parallelization benefit
- Heavy use of statistics and data science
- Visualization is core requirement
- Complex data transformations (finding junctions, aggregations)
- Libraries like `pysam`, `pandas`, `matplotlib` are perfect

**Example Implementation:**

```python
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import json
from collections import defaultdict

# Open BAM file
bam = pysam.AlignmentFile('sample.bam', 'rb')

# Calculate basic statistics
total_reads = 0
mapped_reads = 0
unmapped_reads = 0
quality_scores = []

for read in bam:
    total_reads += 1
    quality_scores.extend(read.query_qualities)
    
    if read.is_unmapped:
        unmapped_reads += 1
    else:
        mapped_reads += 1

# Statistics dictionary
stats = {
    'total_reads': total_reads,
    'mapped_reads': mapped_reads,
    'unmapped_reads': unmapped_reads,
    'mapping_rate': mapped_reads / total_reads * 100,
    'mean_quality': sum(quality_scores) / len(quality_scores),
    'median_quality': sorted(quality_scores)[len(quality_scores)//2]
}

# Find junctions
junctions = defaultdict(int)
junction_list = []

bam.reset()  # Reset to beginning
for read in bam:
    if 'N' in read.cigarstring:
        # This read spans an intron
        junctions[read.reference_name] += 1
        
        # Parse CIGAR to find exact junction coordinates
        cigar = read.cigartuples
        current_pos = read.reference_start
        
        for op, length in cigar:
            if op == 3:  # N in CIGAR = intron
                junction = {
                    'chromosome': read.reference_name,
                    'start': current_pos,
                    'end': current_pos + length,
                    'length': length,
                    'strand': '-' if read.is_reverse else '+'
                }
                junction_list.append(junction)
            elif op in [0, 2]:  # M or D = moves reference position
                current_pos += length

# Coverage analysis
coverage = defaultdict(int)
bam.reset()
for read in bam:
    for pos in range(read.reference_start, read.reference_end):
        coverage[pos] += 1

# Create visualizations
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: Read statistics
ax = axes[0, 0]
categories = ['Mapped', 'Unmapped']
values = [stats['mapped_reads'], stats['unmapped_reads']]
colors = ['green', 'red']
ax.bar(categories, values, color=colors)
ax.set_ylabel('Number of Reads')
ax.set_title(f'Read Mapping ({stats["mapping_rate"]:.1f}% mapped)')
ax.set_yscale('log')

# Plot 2: Quality score distribution
ax = axes[0, 1]
ax.hist(quality_scores, bins=50, edgecolor='black')
ax.set_xlabel('Quality Score')
ax.set_ylabel('Frequency')
ax.set_title(f'Quality Score Distribution (Mean: {stats["mean_quality"]:.1f})')
ax.axvline(stats['mean_quality'], color='red', linestyle='--', label='Mean')
ax.legend()

# Plot 3: Junctions per chromosome
ax = axes[1, 0]
chromosomes = list(junctions.keys())
counts = [junctions[c] for c in chromosomes]
ax.barh(chromosomes, counts)
ax.set_xlabel('Number of Junctions')
ax.set_title('Splicing Junctions by Chromosome')

# Plot 4: Coverage distribution
ax = axes[1, 1]
coverage_values = list(coverage.values())
ax.hist(coverage_values, bins=50, edgecolor='black')
ax.set_xlabel('Coverage Depth')
ax.set_ylabel('Frequency')
ax.set_title('Coverage Distribution')
ax.set_yscale('log')

plt.tight_layout()
plt.savefig('sample_analysis.png', dpi=300)
print("Visualization saved to sample_analysis.png")

# Generate JSON report
report = {
    'sample_statistics': stats,
    'junctions': {
        'total_spliced_reads': len(junction_list),
        'junctions_per_chromosome': dict(junctions),
        'sample_junctions': junction_list[:100]  # First 100 for report
    },
    'coverage': {
        'mean': sum(coverage.values()) / len(coverage),
        'median': sorted(coverage.values())[len(coverage)//2],
        'min': min(coverage.values()),
        'max': max(coverage.values())
    }
}

with open('analysis_report.json', 'w') as f:
    json.dump(report, f, indent=2)

print("Report saved to analysis_report.json")
print(f"Sample statistics: {stats}")
```

**Why NOT Nextflow:**
- Only 1 file (no parallelization benefit)
- Heavy data science work (Python's strength)
- Complex statistical operations
- Visualization (not Nextflow's purpose)
- Would be overkill for single file

**When Python is best:**
- ✅ Single file processing
- ✅ Data science/statistics
- ✅ Complex algorithm development
- ✅ Visualization and plots
- ✅ Scientific computing

**Hybrid approach (if you had many files):**
```groovy
// Nextflow coordinates the work
process ANALYZE_SAMPLE {
    input:
        path bam_file
    output:
        path "*.json"
    script:
        """
        python3 analyze.py ${bam_file}
        """
}

workflow {
    bam_files = Channel.fromPath("*.bam")
    ANALYZE_SAMPLE(bam_files)
}
```

This combines Python (analysis) with Nextflow (parallelization).

---

## Exercise 2: Understanding Parallelization

### Python Sequential Calculation

**Code:**
```python
# Sequential processing - each step must complete before next starts
for sample in samples:  # 1000 samples
    quality_check(sample)     # 5 minutes per sample
    align(sample)             # 30 minutes per sample
    call_variants(sample)     # 15 minutes per sample
```

**Calculation:**
```
Time = 1000 samples × (5 + 30 + 15) minutes
     = 1000 × 50 minutes
     = 50,000 minutes
     = 833 hours
     = 35 days
```

**Timeline visualization:**
```
Sample 1:  QC(5) → ALIGN(30) → VAR(15) ━━━━━━━━━━━━━━━━━━━ 50 min
Sample 2:                        QC(5) → ALIGN(30) → VAR(15) ━━━━ 50 min
Sample 3:                                             QC(5) → ALIGN(30) → VAR(15) ━━ 50 min
...
Sample 1000:                                                                           50 min

Total: 50,000 minutes
```

### Nextflow Parallel Calculation

**Code:**
```groovy
workflow {
    samples = Channel.fromPath("*.fastq")
    
    qc_out = QC(samples)           // All 1000 run in parallel (respecting resources)
    aligned = ALIGN(qc_out)        // Start as soon as each QC completes
    VARIANTS(aligned)              // Start as soon as each ALIGN completes
}
```

**Calculation:**
Assuming you have 8 CPU cores and good I/O, Nextflow can run ~4-8 jobs in parallel:

```
Timeline (8 parallel jobs):
Time 0-5min:    QC for samples 1-8 run in parallel
Time 5-35min:   QC for samples 9-16, ALIGN for samples 1-8 in parallel
Time 35-50min:  QC for samples 17-24, ALIGN for samples 9-16, VAR for samples 1-8 in parallel
...

Total time ≈ (1000 × 50 minutes) / 8 parallel capacity + overhead
          ≈ 50,000 / 8 + ~100 min
          ≈ 6,350 minutes
          ≈ 106 hours
          ≈ 4.4 days
```

**Actual speedup factor:**
```
Python time: 50,000 minutes (35 days)
Nextflow time: 6,350 minutes (4.4 days)
Speedup: 50,000 / 6,350 ≈ 7.9x
```

With a 20-core cluster:
```
Nextflow time: 50,000 / 20 + overhead ≈ 2,600 minutes ≈ 1.8 days
Speedup: 50,000 / 2,600 ≈ 19x
```

With a 100-core cluster:
```
Nextflow time: 50,000 / 100 + overhead ≈ 650 minutes ≈ 0.45 days
Speedup: 50,000 / 650 ≈ 77x
```

**Key insight:** The more resources available, the more Nextflow can parallelize. Python scripts would need rewriting for each environment. Nextflow uses the same code everywhere.

---

## Exercise 3: Recognizing Workflow Patterns

### Original Python Script Analysis

The Python script does:
1. **Find paired-end FASTQ files** - Tumor and normal samples
2. **For each pair:**
   - Call somatic variants with Strelka
   - Annotate variants with VEP
3. **Merge all results** - Combine into one VCF

### Nextflow Workflow Structure

```groovy
// Process 1: Somatic variant calling
process STRELKA {
    input:
        tuple val(sample_id), path(tumor), path(normal)
    output:
        tuple val(sample_id), path("*.vcf")
    script:
        """
        strelka \
            --tumor ${tumor} \
            --normal ${normal} \
            --output ${sample_id}.vcf
        """
}

// Process 2: Variant annotation
process VEP {
    input:
        tuple val(sample_id), path(vcf)
    output:
        tuple val(sample_id), path("*.annotated.vcf")
    script:
        """
        vep \
            --input ${vcf} \
            --output_file ${sample_id}.annotated.vcf \
            --vcf
        """
}

// Process 3: Merge all samples
process MERGE_VARIANTS {
    publishDir "results", mode: 'copy'
    
    input:
        path vcf_files
    output:
        path "merged_variants.vcf"
    script:
        """
        bcftools concat ${vcf_files} -o merged_variants.vcf
        bcftools sort merged_variants.vcf -o merged_variants.sorted.vcf
        mv merged_variants.sorted.vcf merged_variants.vcf
        """
}

// Workflow orchestration
workflow {
    // Create channel of tumor-normal pairs
    sample_pairs = Channel.fromFilePairs("data/*_{tumor,normal}.bam")
    
    // Rename tuple from [sample_id, [tumor_file, normal_file]]
    // to [sample_id, tumor_file, normal_file]
    pairs_reformatted = sample_pairs
        .map { sample_id, files ->
            def tumor = files.find { it.name.contains('tumor') }
            def normal = files.find { it.name.contains('normal') }
            [sample_id, tumor, normal]
        }
    
    // Execute pipeline
    // Strelka runs on all samples in parallel
    variants = STRELKA(pairs_reformatted)
    
    // VEP runs on all variants in parallel
    annotated = VEP(variants)
    
    // Merge waits for all annotation to complete
    MERGE_VARIANTS(annotated.map { id, vcf -> vcf }.collect())
}
```

### Key Differences from Python

| Aspect | Python | Nextflow |
|--------|--------|----------|
| **Explicit loops** | `for tumor, normal in zip()` | Implicit, declarative |
| **Parallelization** | None (sequential) | Automatic |
| **Pair management** | Manual pairing with `zip()` | Using tuples |
| **Error handling** | Manual try-except | Built-in with `errorStrategy` |
| **Resumability** | No restart mechanism | `-resume` flag |
| **Sample tracking** | Lose identity midway | Sample ID stays with data |
| **File management** | Manual output naming | `publishDir` handles it |

### Detailed Workflow Explanation

**Channel operations:**

```groovy
// Input: BAM files like:
// data/sample1_tumor.bam
// data/sample1_normal.bam
// data/sample2_tumor.bam
// data/sample2_normal.bam

// Create paired channel
sample_pairs = Channel.fromFilePairs("data/*_{tumor,normal}.bam")

// This produces:
// [sample1, [tumor.bam, normal.bam]]
// [sample2, [tumor.bam, normal.bam]]

// Reformat for clarity
reformatted = sample_pairs.map { sample_id, files ->
    def tumor = files.find { it.name.contains('tumor') }
    def normal = files.find { it.name.contains('normal') }
    [sample_id, tumor, normal]
}

// Now produces:
// [sample1, sample1_tumor.bam, sample1_normal.bam]
// [sample2, sample2_tumor.bam, sample2_normal.bam]
```

**Parallelization happens automatically:**

```
STRELKA:
Sample1 tumor+normal → Call variants → sample1.vcf
Sample2 tumor+normal → Call variants → sample2.vcf  (parallel)
Sample3 tumor+normal → Call variants → sample3.vcf  (parallel)
...

VEP:
sample1.vcf → Annotate → sample1.annotated.vcf
sample2.vcf → Annotate → sample2.annotated.vcf  (parallel)
sample3.vcf → Annotate → sample3.annotated.vcf  (parallel)
...

MERGE:
(waits for all VEP to complete)
sample1.annotated.vcf + sample2.annotated.vcf + ... → merged_variants.vcf
```

### Why This is Better in Nextflow

1. **Clarity** - The workflow reads like a narrative
2. **Parallelization** - Automatic and optimal
3. **Scalability** - Add 1000 samples, same code, just slower
4. **Fault tolerance** - One sample fails? Restart just that one
5. **Reproducibility** - Same code works on laptop, cluster, cloud
6. **Pair integrity** - Sample IDs stay with their data throughout

---

## Summary: When to Use Each Tool

**Python is best for:**
- ✅ Single file or in-memory data processing
- ✅ Data science and statistical analysis
- ✅ Complex transformations and algorithms
- ✅ Visualization and plotting
- ✅ When you need libraries like pandas, numpy, scipy

**Nextflow is best for:**
- ✅ Coordinating multiple tools
- ✅ Processing many samples/files
- ✅ Need parallelization across resources
- ✅ Building production pipelines
- ✅ Scaling from laptop to cluster/cloud
- ✅ Want built-in resumability

**Use both together:**
```
Nextflow orchestrates → Python processes → External tools execute
(Who does what)        (How things work)   (The actual bioinformatics)
```

---

## Day 1 Complete!

You now understand:
- What Nextflow is (orchestration tool)
- Why it exists (solve parallelization, resumability, portability)
- When to use Python vs Nextflow
- How they work together
- The core concepts (processes, channels, workflows)

You're ready for Day 2: Groovy Essentials!
