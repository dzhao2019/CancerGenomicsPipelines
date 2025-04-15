# Handling multi-mapped reads in RNA-seq

The paper by Deschamps-Francoeur et al. was published in Computational and Structural Biotechnology Journal in 2020, addresses the challenges posed by multi-mapped reads in RNA sequencing (RNA-seq) data analysis. Multi-mapped reads are sequences that align to multiple locations in the genome, complicating accurate quantification of gene and transcript expression levels. This document summarizes the key points from the paper and provides practical recommendations for handling multi-mapped reads in RNA-seq analyses.

## Overview

Eukaryotic genomes often contain duplicated sequences resulting from mechanisms such as recombination, whole-genome duplication, and retrotransposition. These duplications can lead to RNA-seq reads that map to multiple genomic locations, known as multi-mapped reads. Effectively managing these reads is crucial for accurate gene and transcript quantification.

## Mechanisms Leading to Sequence Duplication

### Recombination and Whole-Genome Duplication

Unequal crossing-over during recombination can result in tandem gene duplications.


Whole-genome duplication events can create paralogous gene copies, contributing to sequence redundancy.


### Transposable Elements

Transposons can move within the genome, creating new copies of themselves and leading to sequence duplication.


Retrotransposons, through a "copy and paste" mechanism, can generate processed pseudogenes by reverse transcribing RNA into DNA and inserting it at new locations.


## Impact on RNA-seq Data Analysis

Multi-mapped reads present challenges in accurately quantifying gene and transcript expression levels. The proportion of these reads varies depending on factors such as the organism, sample type, library preparation protocol, and computational pipeline used. Effectively handling multi-mapped reads is essential to avoid biases in downstream analyses.


## Strategies for Handling Multi-Mapped Reads
Several computational approaches have been developed to address the challenges posed by multi-mapped reads:

### Ignoring Multi-Mapped Reads

Excluding reads that map to multiple locations to ensure only uniquely mapped reads are considered.

> Tools: HTSeq-count, STAR geneCounts, featureCounts (default settings).


### Counting Each Alignment

Assigning a count to each alignment of a multi-mapped read, which can inflate expression estimates.

> Tools: featureCounts (with -M option).

### Equal Distribution

Splitting the count of a multi-mapped read equally among all its alignment locations.


> Tools: featureCounts (with -M --fraction options), Cufflinks.

### Rescue Methods

Distributing multi-mapped reads based on the proportion of uniquely mapped reads at each location.


> Tools: Cufflinks (with -u option), CoCo.


### Expectation-Maximization (EM) Algorithms

Using statistical models to probabilistically assign multi-mapped reads based on the overall read distribution.


> Tools: RSEM, IsoEM, EMASE.


## Recommendations

### Assess the Extent of Multi-Mapping

Evaluate the proportion of multi-mapped reads in dataset to inform the choice of handling strategy.

### Choose Appropriate Tools

Select quantification tools that align with analysis goals and can effectively manage multi-mapped reads.

### Document Parameters

Clearly document the parameters and strategies used for handling multi-mapped reads to ensure reproducibility.


## References

Deschamps-Francoeur, G., Simoneau, J., & Scott, M. S. (2020). Handling multi-mapped reads in RNA-seq. Computational and Structural Biotechnology Journal, 18, 1569–1576. 
[Link](https://doi.org/10.1016/j.csbj.2020.06.014)