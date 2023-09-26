## Features:
- Count umi for each gene in each barcode.
- Filter UMI: 
    1. Cannot contain 'N'.
    2. Cannot be a multimer, such as 'AAAAAAAAAA'.
    3. Cannot have base quality lower than 10.


## Arguments
`--probe` Required. Probe sequence file path, 2 columns, gene_name	probe_sequence.

`--fq` Required. Fastq file path.

`--n_mismatch` The maximum allowed mismatch of bases. Default `4`.

`--clean` Remove intermediate files.

