# Tumor RNA-Seq Analysis Using FeatureCounts: Full Workflow from FASTQ to Gene Counts

Before running **FeatureCounts**, you need to process **tumor RNA-Seq FASTQ files** through several essential bioinformatics steps. Below is a detailed step-by-step workflow to prepare your data.

---

## Full Workflow: From FASTQ to FeatureCounts

| **Step** | **Description** | **Tools Used** |
|----------|----------------|---------------|
| **1. Quality Control (QC)** | Check and trim raw RNA-Seq reads | FastQC |
| **2. Read Trimming** | Remove adapters and low-quality bases | Fastp |
| **3. Read Alignment** | Map reads to a reference genome | STAR, Salmon |
| **4. Sorting & Indexing BAM Files** | Prepare BAM for counting | Samtools |
| **5. Run FeatureCounts** | Count mapped reads per gene | FeatureCounts |
| **6. Downstream Analysis** | Normalize, visualize, and analyze gene expression | DESeq2 |

---

## Step 1: Quality Control (QC)

Before running **FeatureCounts**, check the quality of raw **FASTQ** reads.

### **Run FastQC**
```
fastqc tumor_sample_R1.fastq.gz tumor_sample_R2.fastq.gz -o QC_reports/
```

#### Inputs:

* tumor_sample_R1.fastq.gz → Forward reads.
* tumor_sample_R2.fastq.gz → Reverse reads.
* -o QC_reports/ → Save output in QC_reports/ directory.


## Step 2: Read Trimming Using Fastp

Trimming improves alignment accuracy by removing low-quality bases and adapters.

### **Run Fastp**
```
fastp -i tumor_sample_R1.fastq.gz -I tumor_sample_R2.fastq.gz \
      -o trimmed/tumor_sample_R1.trimmed.fastq.gz \
      -O trimmed/tumor_sample_R2.trimmed.fastq.gz \
      --json fastp_report.json --html fastp_report.html
```

#### Inputs:

* tumor_sample_R1.fastq.gz → Forward reads.
* tumor_sample_R2.fastq.gz → Reverse reads.

#### Outputs:
* Trimmed reads in trimmed/
* QC reports: fastp_report.json, fastp_report.html

## Step 3: Read Alignment (Mapping)

Align reads to the human genome (hg38).

### **Using STAR (Recommended)**
```
STAR --runThreadN 8 --genomeDir /path/to/genome_index \
     --readFilesIn trimmed_reads/tumor_sample_R1.fastq.gz trimmed_reads/tumor_sample_R2.fastq.gz \
     --readFilesCommand zcat --outFileNamePrefix mapped/ \
     --outSAMtype BAM SortedByCoordinate
```

#### Inputs:

* trimmed_reads/tumor_sample_R1.fastq.gz → Forward reads.
* trimmed_reads/tumor_sample_R2.fastq.gz → Reverse reads.
* genome_index human genome as reference

#### Outputs:
* Aligned and sorted BAM file (tumor_sample.bam) in zcat format 

## Step 4: Sort and Index BAM File

FeatureCounts requires sorted and indexed BAM files.

### **Sorting  and Index BAM File**
```
samtools sort -@ 8 -o mapped/tumor_sample_sorted.bam mapped/tumor_sample.bam

samtools index mapped/tumor_sample_sorted.bam
```

#### Outputs:
* tumor_sample_sorted.bam
* tumor_sample_sorted.bam.bai

## Step 5: Count Reads Using FeatureCounts

Use FeatureCounts to assign reads to genes.

### **Run FeatureCounts**
```
featureCounts -T 8 -a /path/to/annotation.gtf   \
              -o counts.txt mapped/tumor_sample_sorted.bam
              
featureCounts -p -T 8 -s 1  \
              -a /path/to/annotation.gtf   \
              -o  sample1.featureCounts.txt sample1.bam sample2.bam
              
```

#### Options Explained:
* -p → If paired-end sequencing.
* -T 8 → Uses 8 threads for faster processing.
* -a annotation.gtf → Annotation file (e.g., GENCODE, Ensembl).
* -o counts.txt → Output gene count matrix.
* -s $strandedness → Defines the library strandedness:
  - 0 → Unstranded
  - 1 → Forward stranded
  - 2 → Reverse stranded

#### Inputs:
* All BAM files as input.
  - sample1.bam 
  - sample2.bam

#### Output (counts.txt):

```
Geneid   Chr   Start   End   Strand   Length   Tumor_Sample
Gene1    chr1  1000    2000  +        1000      50
Gene2    chr2  3000    4000  -        1000      30

```

## Step 6: Differential Expression Analysis

FeatureCounts output is used for DESeq2.

### **Load Data into DESeq2 (R)**
```
library(DESeq2)
counts <- read.table("counts.txt", header=T, row.names=1)
colData <- data.frame(row.names = colnames(counts), condition = c("tumor", "normal"))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
```


