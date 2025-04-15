# **RNA-Seq Data Analysis for Cancer Sequencing**

This is a structured workflow for analyzing **RNA-Seq data from cancer samples** using **LSF job scheduling**, **FeatureCounts**, and other bioinformatics tools. The pipeline ensures efficient processing, alignment, quantification, and interpretation of cancer RNA-Seq datasets.

---

## **Table of Contents**
1. [Introduction](#introduction)
2. [Workflow Overview](#workflow-overview)
3. [Step-by-Step Guide](#step-by-step-guide)
   - [Step 1: Data Preprocessing](#step-1-data-preprocessing)
   - [Step 2: Alignment with STAR](#step-2-alignment-with-star)
   - [Step 3: BAM Processing with Samtools](#step-3-bam-processing-with-samtools)
   - [Step 4: Read Quantification with FeatureCounts](#step-4-read-quantification-with-featurecounts)
   - [Step 5: Downstream Analysis](#step-5-downstream-analysis)
4. [LSF Job Management](#lsf-job-management)
5. [References](#references)

---

## **Introduction**
RNA-Seq is a powerful technique for studying **gene expression profiles in cancer**. This workflow processes RNA sequencing data from **raw reads to gene expression quantification**, using tools optimized for **high-performance computing (HPC) clusters**.

---

## **Workflow Overview**
The pipeline follows these steps:

1. **Preprocessing**: Quality control and trimming of raw FASTQ reads.
2. **Alignment**: Mapping reads to a reference genome using **STAR**.
3. **BAM Processing**: Sorting, indexing, and filtering aligned reads using **Samtools**.
4. **Quantification**: Counting gene-level expression using **FeatureCounts**.
5. **Downstream Analysis**: Normalization, differential expression, and visualization.

---

## **Step-by-Step Guide**
### **Step 1: Data Preprocessing**
- **Quality control**: Assess raw reads using **FastQC**.
- **Trimming**: Remove adapters and low-quality bases with **Fastp**.

**Example LSF Job Submission:**
```
bsub -J fastqc -o fastqc.out -e fastqc.err -n 4 -R "rusage[mem=8000]" "fastqc -t 4 sample_R1.fastq.gz sample_R2.fastq.gz"
```
### **Step 2: Alignment with STAR**

* Genome index setup (if not already built).
* Align reads to the reference genome.

Example LSF Job Submission:
```
bsub -J star_align -o star.out -e star.err -n 8 -R "rusage[mem=32000]" -W 6:00 \
"STAR --genomeDir /path/to/star_index --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
--runThreadN 8 --readFilesCommand zcat --outFileNamePrefix output/"

```

### **Step 3: BAM Processing with Samtools**
Sort and index BAM files for downstream analysis.

Example LSF Job Submission:
```
bsub -J samtools -o samtools.out -e samtools.err -n 4 -R "rusage[mem=16000]" -W 2:00 \
"samtools sort -@ 4 -o output.sorted.bam output.Aligned.unsort.bam && samtools index output.sorted.bam"

```

### **Step 4: Read Quantification with FeatureCounts**
Count mapped reads per gene using FeatureCounts.
See detailed [FeatureCounts Guide](https://github.com/dzhao2019/BioPipelineCraft/blob/main/RNA/FeatureCount.md).

Example LSF Job Submission:
```
bsub -J featureCounts -o featureCounts.out -e featureCounts.err -n 8 -R "rusage[mem=16000]" -W 4:00 \
"featureCounts -T 8 -a annotation.gtf -o counts.txt output.sorted.bam"
```


### **Step 5: Downstream Analysis**

* Differential Expression Analysis (DESeq2, edgeR).
* Visualization (PCA, heatmaps, volcano plots).

## **LSF Job Management**
For job submission and monitoring on LSF clusters, see the [LSF Command Summary](https://github.com/dzhao2019/BioPipelineCraft/blob/main/RNA/LSF%20Command%20Summary.md).

Essential Commands:

```
bjobs             # Check job status
bkill <jobID>     # Cancel a job
bqueues -l        # View available queues
```

## **References**

- [STAR: Spliced Transcripts Alignment to a Reference](https://github.com/alexdobin/STAR)  
- [FeatureCounts: Read Quantification](http://subread.sourceforge.net/)  
- [FastQC: Quality Control](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
- [LSF User Guide](https://www.ibm.com/docs/en/spectrum-lsf)  
