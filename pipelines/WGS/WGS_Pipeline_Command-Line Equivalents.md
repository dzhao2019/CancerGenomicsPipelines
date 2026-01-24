# WGS Nextflow Pipeline Command-Line Equivalents

These commands are derived from the script blocks within each process, stripped of Nextflow-specific syntax (e.g., variables like ${params.ref}), and presented as standalone Bash commands that can be run in a shell to verify the tools’ execution. Assume typical input/output file names based on the pipeline’s logic and include module loads where specified, so the environment can be replicated. 

________________________________________

## Assumptions

<p>• Environment: Commands assume a Linux system with tools installed or loaded via modules (e.g., module load as in WGS script).
<p>• Parameters: params values are substituted with their defaults in nextflow script (e.g., params.ref = "ref_hg38/UCSC_hg38.fa").
<p>• Inputs: using example file names reflecting the pipeline’s flow (e.g., sample1_R1.merged.fastq.gz).

________________________________________


## Pipeline and Command-Line Equivalents

### Workflow 1: Data Pre-processing
#### i. Merge lanes 

Task: Merges FASTQ lane files into R1 and R2.

```commandline
cat sample1_lane1_R1.fastq.gz sample1_lane2_R1.fastq.gz > sample1_merged_R1.fastq.gz
cat sample1_lane1_R2.fastq.gz sample1_lane2_R2.fastq.gz > sample1_merged_R2.fastq.gz
```
#### ii. Sample Quality control (QC)

**FastQC**: Quality control tool for sequencing reads.

```commandline
# Initial FastQC
module load openjdk/11
module load fastqc/0.12.1

# Run FastQC
fastqc sample1_merged_R1.fastq.gz sample1_merged_R2.fastq.gz
```


#### iii.Quality Trimming: Adapter and quality trimming

**Fastp**: Trims adapters (auto-detected) and low-quality bases.

```commandline
# Run Fastp
# Trimmed FASTQs, HTML/JSON reports.
module load fastp/0.19.5
fastp -i sample1_merged_R1.fastq.gz -I sample1_merged_R2.fastq.gz \
      -o sample1.trimmed.R1.fastq.gz -O sample1.trimmed.R2.fastq.gz \
      -l 30 -W 32 --html sample1.fastp.html --json sample1.fastp.json
```


**Parameters for Fastp:**

* -i: Input R1 FASTQ file (sample1_merged_R1.fastq.gz).
* -I: Input R2 FASTQ file (sample1_merged_R2.fastq.gz).
* -l 30: Minimum read length after trimming (30 bp).
  - Purpose: Discards reads shorter than 30 bp post-trimming to ensure quality.
* -W 32: Sliding window size for quality trimming (32 bp).
  - Purpose: Checks quality in 32-bp windows; trims if average quality drops below a threshold (default Q20 unless overridden).
* -o: Output trimmed R1 (sample1.trimmed.R1.fastq.gz).
* -O: Output trimmed R2 (sample1.trimmed.R2.fastq.gz).
* --html sample1.fastp.html: HTML report of trimming stats.
* --json sample1.fastp.json: JSON report for programmatic use.

**Interpretation:**

Trims adapters (auto-detected) and low-quality bases, keeping reads ≥30 bp. Window size of 32 balances sensitivity and speed.

#### iv. Fastqc on trimmed fastq read

```commandline
# Run FastQC on Trimmed fastq     
module load fastqc/0.12.1
fastqc sample1.trimmed.R1.fastq.gz sample1.trimmed.R2.fastq.gz     
```

#### v. MultiQC on pre-processed sample

**Task**: Runs MultiQC on merged and trimmed FastQC/Fastp outputs.

1. MultiQC on Merged

```commandline
module load multiqc/1.19
multiqc sample1_merged_R1_fastqc.zip sample1_merged_R2_fastqc.zip
mv multiqc_data merged_multiqc_data
mv multiqc_report.html merged_multiqc.html
```
<p>Default output: multiqc_report.html
<p>Renames outputs for clarity (e.g., merged_multiqc.html)


2. MultiQC on Trimmed:
```commandline
# Includes Fastp JSON alongside FastQC ZIPs.
multiqc sample1.trimmed.R1_fastqc.zip sample1.trimmed.R2_fastqc.zip sample1.fastp.json
mv multiqc_data trimmed_multiqc_data
mv multiqc_report.html trimmed_multiqc.html
```
<p>Renames outputs for clarity (e.g., trimmed_multiqc.html)

### Workflow 2: Genome Alignment with BWA mem

**Task**: Aligns trimmed reads to the genome with BWA-MEM, sorts with Samtools, and runs Qualimap.


  Input: trimmed.fastq.gz
  Output: sorted bam

```commandline
/path/to/bwa mem ref_hg38/UCSC_hg38.fa sample1.trimmed.R1.fastq.gz sample1.trimmed.R2.fastq.gz \
             -t 20 -R "@RG\tID:sample1\tPL:Illumina\tPU:None\tLB:None\tSM:sample1" \
              | samtools sort -o sample1.sorted.bam
samtools index sample1.sorted.bam
```

Parameters:
* bwa mem: Aligns reads to a reference.
* ref_hg38/UCSC_hg38.fa: Reference genome FASTA.
* sample1.trimmed.R1.fastq.gz sample1.trimmed.R2.fastq.gz: Input paired-end reads.
* -t 20: Number of threads (20).
  * Purpose: Speeds up alignment using 20 CPU cores.
* -R "@RG\tID:sample1\tPL:Illumina\tPU:None\tLB:None\tSM:sample1": Read group header.
  * Purpose: Adds metadata (ID, platform, sample name) to BAM for downstream tools like GATK.
* | samtools sort: Pipes BWA output (SAM) to sort into BAM.
* -o sample1.sorted.bam: Output sorted BAM.

**Interpretation**:
<p>Aligns reads to hg38, adds read group info, and sorts by coordinate for efficiency.
<p>Creates sample1.sorted.bam.bai index for fast BAM access.

### Workflow 3: Mapping statistics with Qualimap

**Qualimap** analyzes BAM files and generates comprehensive quality control (QC) metrics. 
Its primary mode, **BAMQC** focuses on assessing how well sequencing reads align to a reference genome. It answers questions like:
<p>Are the reads evenly distributed across the genome?
<p>How much of the genome is covered?
<p>Are there biases or artifacts affecting the data?

```commandline
module load openjdk/11
module load qualimap/2.3
qualimap bamqc --java-mem-size=10G sample1.sorted.bam \
               -c -nw 400 -hm 3 -outdir sample1_raw_qualimap \
               -outfile sample1_raw.pdf -outformat PDF
```

**Parameters**:
* bamqc: BAM quality control.
* --java-mem-size=10G: Allocates 10 GB RAM.
  * Purpose: Ensures Qualimap handles large BAMs.
* sample1.sorted.bam: Input BAM.
* -c: Counts reads per region.
  * Purpose: Includes coverage stats.
* -nw 400: Number of windows (400).
  * Purpose: Divides genome into 400 bins for coverage analysis.
* -hm 3: Histogram bins (3).
  * Purpose: Simplifies coverage histogram.
* -outdir sample1_raw_qualimap: Output directory.
* -outfile sample1_raw.pdf: Output file name.
* -outformat PDF: PDF report.

**Interpretation**: Assesses alignment quality (coverage, mapping stats).

### Workflow 4: Pre-processing of sorted bam file with gatk

#### i. Markduplicates reads.
**Task**: Pre-processes BAMs with GATK and runs MarkDuplicates

Identifies PCR duplicates to improve variant calling accuracy.

```commandline
module load openjdk/17
module load gatk/4.5.0.0.8
mkdir gatk_temp
gatk --java-options "-Xmx32G -Djava.io.tmpdir=gatk_temp" MarkDuplicates \
     -I sample1.sorted.bam -O sample1.dup.bam -M sample1.metrics.txt
```

**Parameters**:
* --java-options "-Xmx32G -Djava.io.tmpdir=gatk_temp": JVM settings.
* -Xmx32G: 32 GB RAM.
* -Djava.io.tmpdir=gatk_temp: Temp directory.
* MarkDuplicates: Marks duplicate reads.
* -I sample1.sorted.bam: Input BAM.
* -O sample1.dup.bam: Output BAM with duplicates marked.
* -M sample1.metrics.txt: Metrics file.

#### ii. SetNmMdAndUqTags

**Task**: Corrects NM (edit distance), MD (mismatch), and UQ (quality) tags for downstream GATK tools

```commandline
gatk --java-options "-Xmx32G -Djava.io.tmpdir=gatk_temp" SetNmMdAndUqTags \
     -R ref_hg38/UCSC_hg38.fa -I sample1.dup.bam -O sample1.markdup.bam
```

**Parameters**:
* -R ref_hg38/UCSC_hg38.fa: Reference genome.
* -I sample1.dup.bam: Input BAM.
* -O sample1.markdup.bam: Output BAM with fixed tags.


#### iii. Baserecaliberation (BQSR)
**Task**: Calibrates base quality scores using known variants to reduce sequencing errors.

```commandline
gatk --java-options "-Xmx32G -Djava.io.tmpdir=gatk_temp" BaseRecalibrator \
     --input sample1.markdup.bam -R ref_hg38/UCSC_hg38.fa \
     --known-sites 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
     --known-sites Homo_sapiens_assembly38.dbsnp138.vcf \
     --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
     -O sample1.recal.table
```
**Parameters**:
* --input sample1.markdup.bam: Input BAM.
* -R ref_hg38/UCSC_hg38.fa: Reference genome.
* --known-sites: Known variant databases.
* 1000G_phase1.snps...: SNPs from 1000 Genomes.
* Homo_sapiens_assembly38.dbsnp138.vcf: dbSNP variants.
* Mills_and_1000G_gold_standard.indels...: Indels from Mills/1000G.
* -O sample1.recal.table: Recalibration table.

#### iv. Apply BQSR
**Task**: Applies BQSR adjustments to improve variant calling.

```commandline
gatk --java-options "-Xmx32G -Djava.io.tmpdir=gatk_temp" ApplyBQSR \
     -R ref_hg38/UCSC_hg38.fa -I sample1.markdup.bam -O sample1.ready.bam \
      -bqsr sample1.recal.table
```
**Parameters**:
* -R ref_hg38/UCSC_hg38.fa: Reference genome.
* -I sample1.markdup.bam: Input BAM.
* -O sample1.ready.bam: Output recalibrated BAM.
* -bqsr sample1.recal.table: Recalibration table.

#### v. Deduplicated statistics with Qualimap
**Task**: Evaluates BAM post-duplicate marking.
Runs the BAM QC mode of Qualimap to evaluate alignment quality


```commandline
module load openjdk/11
module load qualimap/2.3
#Evaluates BAM post-duplicate marking
qualimap bamqc --java-mem-size=10G sample1.markdup.bam \
               -sd -c -nw 400 -hm 3 \
               -outdir sample1_markdup_qualimap \
               -outfile sample1_markdup_qualimap.pdf -outformat PDF
               
#Assesses final BAM quality after BQSR.
qualimap bamqc --java-mem-size=10G sample1.ready.bam \
               -sd -c -nw 400 -hm 3 \
               -outdir sample1_bam_ready_qualimap \
               -outfile sample1_bam_ready_qualimap.pdf -outformat PDF               
```

**Parameters**:
* --java-mem-size=10G: Allocates 10 GB of RAM to the Java process.
* -sd: (skip duplicates) to exclude marked duplicates from stats.
  * Purpose: Skips duplicate reads (flagged with 0x400 in the BAM) in the analysis. Focuses stats on unique reads only, reflecting the effective data after deduplication. 
* -c: Counts reads per region (genome-wide unless a GFF is provided).
* -nw 400: Sets the number of windows to 400 for coverage analysis.
  * Purpose: Divides the genome into 400 bins to compute coverage distribution. A moderate number balances granularity and computation time.
* -hm 3: Uses 3 bins for the coverage histogram
  * Purpose: Simplifies the coverage plot into three categories (e.g., low, medium, high), making it easier to interpret.


### Workflow 5: Mutect2 variants calling.

**Task**: Calls variants with GATK Mutect2.
Calls somatic variants by comparing tumor and normal samples

```commandline
module load bioinfo/gatk/4.1.8
module load java/openjdk/15
gatk Mutect2 -R ref_hg38/UCSC_hg38.fa \
            -I sample1.ready.bam -I sample1_normal.ready.bam \
            -normal sample1_normal -O sample1.vcf.gz
```
**Parameters**:
* -R ref_hg38/UCSC_hg38.fa: Reference genome.
* -I sample1.ready.bam: Tumor BAM.
* -I sample1_normal.ready.bam: Normal BAM.
* -normal sample1_normal: Normal sample ID (matches read group SM).
  * Purpose: Distinguishes tumor vs. normal for somatic calling.
* -O sample1.vcf.gz: Output VCF (gzipped).

