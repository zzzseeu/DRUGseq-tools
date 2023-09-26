## Features:
- Count umi for each gene in each barcode.
- Filter UMI: 
    1. Cannot contain 'N'.
    2. Cannot be a multimer, such as 'AAAAAAAAAA'.
    3. Cannot have base quality lower than 10.

## Outputs:
- `{sample}.detail.txt`: UMI, read count raw file. 4 columns, Barcode       geneID  UMI     count.
- `{sample}.features.txt`: Gene features record file. 3 columns, gene_name  gene_id gene_biotype.
- `{samples}.metadata.txt`: Sample meta data. 5 columns, Barcode    readcount       UMI2(UMI counts>2)      UMI     GeneCount.
- `{sample}.matrix.txt`: Gene expression matrix. Columns are genes, rows are samples.
- `{sample}.matrix.h5ad`: Gene expression matrix h5 file for downstream analysis.


## Arguments
`--bam` Sorted featureCounts output bamfile.

`--gtf` GTF file path.

`--outTable` Output gene expression count table.

