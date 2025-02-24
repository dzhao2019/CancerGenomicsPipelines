## **Comprehensive Guide to Cancer RNA-Seq Data Analysis Pipeline**

This guide provides a step-by-step explanation of the tools used in the RNA-Seq data analysis pipeline for cancer research. Each tool's purpose, function, and example command lines are described in the order of operation.

---

### **1. Subsampling Raw Sequencing Data**

#### **Tool: FQ_SUBSAMPLE**
- **Purpose:** Randomly subsample FASTQ files to create smaller datasets.
- **Function:** Useful for quick quality control checks, benchmarking pipelines, and normalizing read depths.

#### **Commands:**
- **Subsample by Read Count:**
  ```
  fq_subsample -i sample_R1.fastq -o subsampled_R1.fastq -n 1000000
  ```
- **Subsample by Proportion:**
  ```
  fq_subsample -i sample_R1.fastq -o subsampled_R1.fastq -p 0.1
  ```
- **Output:** Subsampled FASTQ files for downstream QC or testing.

---

### **2. Quality Control of Raw Sequencing Data**

#### **Tool: FastQC**
- **Purpose:** Assess the quality of raw sequencing reads.
- **Function:** Provides metrics like per base quality, GC content, and adapter contamination.

#### **Command:**
```
fastqc sample_R1.fastq sample_R2.fastq -o fastqc_results/
```
- **Output:** HTML and text reports summarizing quality metrics.

---

#### **Tool: fastp**
- **Purpose:** Perform quality control and preprocessing of raw FASTQ files.
- **Function:** Integrates quality filtering, trimming, and adapter removal in one step.

#### **Command:**
```
fastp -i sample_R1.fastq -I sample_R2.fastq -o trimmed_R1.fastq -O trimmed_R2.fastq \
      -h fastp_report.html -j fastp_report.json
```
- **Output:** Trimmed FASTQ files with an HTML report summarizing quality metrics.

#### **Example log file location:**
```
/research/project/TerryFoxRI/G99_RNA1/fastp/log/sample_1_A.fastp.log
```

#### **Example execution with parameters:**
```
fastp --in1 sample_1_A_1.fastq.gz --in2 sample_1_A_2.fastq.gz \
      --out1 sample_1_A_1.fastp.fastq.gz --out2 sample_1_A_2.fastp.fastq.gz \
      --json sample_1_A.fastp.json --html sample_1_A.fastp.html \
      --thread 16 --detect_adapter_for_pe
```

---

## **Post-Trimming QC (FastQC Again)**

### **Purpose**
- Verify improvements in data quality after trimming.

### **Command**
```
fastqc trimmed_R1.fastq trimmed_R2.fastq -o fastqc_trimmed_results/
```

---

## **3. Read Alignment to the Reference Genome**

### **Data Transformation from Trimmed FASTQ to FASTA**

#### **Purpose**
- Convert trimmed FASTQ files to FASTA format for assembly or other sequence-based analyses.

#### **Tools**
- **seqtk**: A lightweight tool for FASTQ to FASTA conversion.

#### **Command**
```
seqtk seq -a trimmed_R1.fastq > trimmed_R1.fasta
seqtk seq -a trimmed_R2.fastq > trimmed_R2.fasta
```

#### **Output**
- FASTA files: `trimmed_R1.fasta`, `trimmed_R2.fasta`

---

### **Assembly Tools for FASTA Files**

#### **Purpose**
- Assemble short reads into longer contiguous sequences (contigs).

#### **Common Assembly Tools**
- **Trinity**: Designed for RNA-Seq data to assemble transcriptomes.

#### **Command**
```
Trinity --seqType fa --left trimmed_R1.fasta --right trimmed_R2.fasta --CPU 8 --max_memory 50G --output trinity_out_dir
```

- **SPAdes**: Versatile assembler for both genomic and transcriptomic data.

#### **Command**
```
spades.py -1 trimmed_R1.fasta -2 trimmed_R2.fasta -o spades_output
```

#### **Output**
- Assembled contigs in FASTA format.

---

### **Annotation Tools for FASTA/FA Files**

#### **Purpose**
- Annotate assembled sequences to identify genes, features, and functions.

#### **Common Annotation Tools**
- **Prokka**: Rapid prokaryotic genome annotation tool.

#### **Command**
```
prokka assembled_contigs.fasta --outdir prokka_output --prefix sample
```

- **BLAST (Basic Local Alignment Search Tool)**: For sequence similarity searches against known databases.

#### **Command**
```
blastn -query assembled_contigs.fasta -db nt -out results.out -evalue 1e-5 -outfmt 6
```

- **InterProScan**: For functional analysis and classification of proteins.

#### **Command**
```
interproscan.sh -i proteins.fasta -o interpro_output.tsv -f tsv
```

#### **Output**
- Annotated files with gene functions, sequence features, and similarity matches.

---

## **Process: PREPROCESS_TRANSCRIPTS_FASTA_GENCODE**

### **Software**
- `sed`

### **Purpose**
- Preprocess transcript FASTA files from GENCODE to ensure proper formatting for downstream analyses.

### **Function**
- Modifies headers, removes unnecessary information, or adjusts formatting in FASTA files to meet pipeline requirements.

### **Command**
```
sed 's/|.*//g' gencode_transcripts.fa > processed_transcripts.fa
```

### **Output**
- Processed FASTA file: `processed_transcripts.fa` (ready for alignment and quantification).

---

## **Tool: STAR (Spliced Transcripts Alignment to a Reference)**

### **Purpose**
- Align trimmed reads to the reference genome.

### **Function**
- Identifies exon-exon junctions in RNA-Seq data.

### **Reference Indexing**
```
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir genome_index \
     --genomeFastaFiles reference_genome.fa --sjdbGTFfile annotation.gtf --sjdbOverhang 100
```

### **Alignment Command**
```
STAR --runThreadN 8 --genomeDir genome_index \
     --readFilesIn trimmed_R1.fastq trimmed_R2.fastq \
     --outFileNamePrefix aligned_output/
```

### **Output**
- Aligned BAM files.

---

## **4. Post-Alignment Processing**

### **Tool: GTF2BED**

#### **Purpose**
- Convert gene annotation files from GTF to BED format.

#### **Function**
- Enables compatibility with tools like BEDTools for coverage analysis and genomic region processing.

#### **Command**
```
gtf2bed < annotation.gtf > annotation.bed
```

#### **Output**
- BED file: `annotation.bed` (for downstream coverage analysis).

---

### **Tool: Picard MarkDuplicates**

#### **Purpose**
- Identify and mark duplicate reads in BAM files.

#### **Function**
- Removes PCR duplicates, which can introduce bias in downstream analyses like variant calling or gene expression quantification.

#### **Command**
```
picard MarkDuplicates I=sorted_aligned.bam O=dedup_aligned.bam M=dedup_metrics.txt
```

#### **Output**
- Deduplicated BAM file: `dedup_aligned.bam`
- Metrics file: `dedup_metrics.txt` (summarizing duplication statistics).

#### **Example Log Location**
```
/research/project/TerryFoxRI/G99_RNA1/star_salmon/picard_metrics/sample_1_A.markdup.sorted.MarkDuplicates.metrics.txt
```

---

## **Tool: Samtools**

### **SAMTOOLS_SORT**
#### **Purpose**
- Sort BAM files by genomic coordinates.

#### **Command**
```
samtools sort -o sorted_aligned.bam aligned_output/Aligned.out.sam
```

#### **Output**
- Sorted BAM file: `sorted_aligned.bam`

---

### **SAMTOOLS_INDEX**
#### **Purpose**
- Index sorted BAM files for efficient data retrieval.

#### **Command**
```
samtools index sorted_aligned.bam
```

#### **Output**
- BAM index file: `sorted_aligned.bam.bai`

---

### **SAMTOOLS_FLAGSTAT**
#### **Purpose**
- Provide statistics on alignment quality and read mapping.

#### **Command**
```
samtools flagstat sorted_aligned.bam
```

#### **Output**
- Summary statistics of the BAM file, including total reads, mapped reads, and duplicates.

---

### **SAMTOOLS_IDXSTATS**
#### **Purpose**
- Display alignment statistics per reference sequence.

#### **Command**
```
samtools idxstats sorted_aligned.bam
```

#### **Output**
- Statistics for each chromosome, including mapped and unmapped read counts.

---

### **SAMTOOLS_STATS**
#### **Purpose**
- Generate detailed statistics for alignment quality.

#### **Command**
```
samtools stats sorted_aligned.bam > alignment_stats.txt
```

#### **Output**
- Comprehensive alignment statistics in `alignment_stats.txt`.

---

## **Tool: BEDTools Genomecov**

#### **Purpose**
- Calculate genome-wide or region-specific coverage from BAM files.

#### **Function**
- Determines the depth of coverage at each base or over specified regions, helping identify areas with high or low read coverage.

#### **Command**
```
bedtools genomecov -ibam sorted_aligned.bam -bg > coverage.bedgraph
```

#### **Output**
- BEDGraph file: `coverage.bedgraph` (showing the coverage depth across the genome).

---

## **Tool: UCSC_BEDCLIP**

#### **Purpose**
- Clip BED intervals that extend beyond the chromosome sizes.

#### **Function**
- Ensures that BED intervals do not exceed the reference genome boundaries.

#### **Command**
```
bedClip coverage.bedgraph chrom.sizes clipped_coverage.bedgraph
```

#### **Output**
- Clipped BEDGraph file: `clipped_coverage.bedgraph` (ensuring valid genomic coordinates).

---

## **Tool: UCSC_BEDGRAPHTOBIGWIG**

#### **Purpose**
- Convert BEDGraph files to BigWig format for efficient visualization.

#### **Function**
- Transforms BEDGraph files into BigWig format, suitable for genome browsers.

#### **Command**
```
bedGraphToBigWig clipped_coverage.bedgraph chrom.sizes coverage.bw
```

#### **Output**
- BigWig file: `coverage.bw` (for visualization in genome browsers).

---

## **Additional Coverage Analysis**

### **Tool: BEDTools Genomecov**
#### **Purpose**
- Calculate genome-wide or region-specific coverage from BAM files.

#### **Function**
- Determines the depth of coverage at each base or over specified regions, helping identify high or low read coverage.

#### **Command**
```
bedtools genomecov -ibam sorted_aligned.bam -bg > coverage.bedgraph
```

#### **Output**
- BEDGraph file: `coverage.bedgraph` (showing coverage depth across the genome).

---

### **Gene Body Coverage Analysis**
#### **Command**
```
geneBody_coverage.py -r annotation.bed -i sorted_aligned.bam -o coverage_output
```

#### **Output**
- Coverage statistics in `coverage_output`.

---

## **5. Gene Expression Quantification**

### **Tool: SALMON_QUANT (salmon)**

#### **Purpose**
- Quantify transcript-level expression from RNA-Seq data.

#### **Function**
- Uses quasi-mapping for fast and accurate transcript quantification.

#### **Command**
```
salmon quant -i transcript_index -l A -1 trimmed_R1.fastq -2 trimmed_R2.fastq -p 8 -o salmon_output
```

#### **Output**
- Quantification files with transcript-level counts and TPM values.

---

### **Tool: SALMON_SE_GENE (bioconductor-summarizedexperiment)**

#### **Purpose**
- Organize gene-level expression data in a structured format.

#### **Function**
- Summarizes transcript-level estimates to gene-level using R.

#### **R Script Example**
```r
library(SummarizedExperiment)
se <- SummarizedExperiment(assays=list(counts=counts, tpm=tpm))
saveRDS(se, file="gene_level_data.rds")
```

#### **Output**
- RDS file with summarized gene-level expression data.

---

### **Tool: SALMON_TX2GENE (python)**

#### **Purpose**
- Map transcript IDs to gene IDs for downstream analysis.

#### **Function**
- Converts transcript-level quantifications to gene-level.

#### **Python Script Example**
```python
import pandas as pd

tx2gene = pd.read_csv('tx2gene_mapping.csv')
salmon_data = pd.read_csv('salmon_output/quant.sf', sep='\t')
merged_data = pd.merge(salmon_data, tx2gene, left_on='Name', right_on='transcript_id')
merged_data.to_csv('gene_counts.csv', index=False)
```

#### **Output**
- CSV file mapping transcripts to genes.

---

### **Tool: SALMON_TXIMPORT (bioconductor-tximeta, r-base)**

#### **Purpose**
- Import and summarize Salmon quantification data for differential expression analysis.

#### **Function**
- Aggregates transcript-level estimates to gene-level using metadata.

#### **R Script Example**
```r
library(tximeta)
files <- list.files("salmon_output", full.names = TRUE)
coldata <- data.frame(files=files, names=c("Sample1", "Sample2"))
txi <- tximeta(coldata)
write.csv(assay(txi, "counts"), file="tximport_gene_counts.csv")
```

#### **Output**
- CSV file with aggregated gene-level counts for DE analysis.

---

### **Tool: FeatureCounts**

#### **Purpose**
- Count the number of reads mapped to each gene.

#### **Function**
- Generates a count matrix for downstream differential expression analysis.

#### **Command**
```
featureCounts -T 8 -a annotation.gtf -o counts.txt sorted_aligned.bam
```

#### **Output**
- A text file with raw counts of reads per gene.

---

## **6. Differential Expression Analysis**

### **Tool: DESeq2 (R Package)**

#### **Purpose**
- Identify differentially expressed genes between conditions.

#### **Function**
- Normalizes count data and performs statistical testing.

#### **R Script Example**
```r
library(DESeq2)
countData <- read.table("counts.txt", header=TRUE, row.names=1)
colData <- data.frame(condition=c("control", "treated"))
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
write.csv(as.data.frame(res), file="differential_expression_results.csv")
```

---

## **7. Functional Enrichment Analysis**

### **Tool: ClusterProfiler (R Package)**

#### **Purpose**
- Perform Gene Ontology (GO) and pathway enrichment analysis.

#### **R Script Example**
```r
library(clusterProfiler)
gene_list <- res$log2FoldChange
names(gene_list) <- rownames(res)
enriched_go <- enrichGO(gene=names(gene_list), OrgDb=org.Hs.eg.db, keyType="ENSEMBL", ont="BP", pAdjustMethod="BH")
dotplot(enriched_go)
```

---

### **Tool: GSEA (Gene Set Enrichment Analysis)**

#### **Purpose**
- Identify pathways significantly enriched in differentially expressed genes.

#### **Command (CLI)**
```
gsea-cli.sh -res expression_data.gct -cls phenotype.cls -gmx gene_sets.gmt -out gsea_results/
```

---

## **8. Data Visualization and Reporting**

### **Tool: MultiQC**

#### **Purpose**
- Aggregate QC results from all analysis steps into a single report.

#### **Command**
```
multiqc . -o multiqc_report/
```

#### **Output**
- An interactive HTML report summarizing the results from all tools.

---

This document continues with further steps such as visualization, reporting, and validation.