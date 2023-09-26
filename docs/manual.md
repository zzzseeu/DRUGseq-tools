## Introduction
DRUGseq-tools is a collection of bioinfomatics analysis pipelines developed at XinRui EMED to process large-scale transcriptome sequencing data generated with Plastech-seq and PHD-seq. These pipelines take paired-end FASTQ files as input and generate output files which can be used for downstream data analysis as well as a summary of QC criteria.

Currently, DRUGseq-tools includes the follwing pipelines:

- `drug plastech` for Molecular tag-based mRNA sequencing data. It performs preprocessing, genome alignment, feature counting, expression matrix generation.

- `drug phd` for targeted transcriptome sequencing data based on molecular tags. It performs preprocessing, UMI filtering and UMI counting.

## [Quick start](quick_start.md)

## Plastech-seq: Transcriptome sequencing analysis based on barcode and UMI.
- [barcode](plastech/barcode.md)
- [trim](plastech/trim.md)
- [mapping](plastech/mapping.md)
- [featureCounts](plastech/featureCounts.md)
- [count](plastech/count.md)
- [analysis](toolkits/analysis.md)
## PHD-seq: Target transcriptome sequencing analysis based on barcode and UMI.
- [preprocess](phd/preprocess.md)
- [count](phd/count.md)
- [analysis](toolkits/analysis.md)
