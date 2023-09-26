## Features: 
- Extract the Barcode and UMI information in R1, and put them in the header of R2 read.
- Filter barcode: Only one base mismatch is allowed at most, 
                  and the mismatched base must be a low-quality base.
                  
## Outputs:
- `{sample}.fq`: R2 fastq data after demultiplexing with barcode and UMI in the seq name. 
    Multi samples will only output one fastq file. 
- `{sample}_summary.txt`: Barcode step summary file.


## Arguments
`--plate_tsv` Required. A record file, 3 columns, seperated by tab:         Library	Path	Info.

`--barcode_range` Barcode position in Read 1. Default `1, 10`.

`--umi_range` UMI position in Read 1. Default `11, 20`.

