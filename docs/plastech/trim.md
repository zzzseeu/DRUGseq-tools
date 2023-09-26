## Features:
- Trim the adapter sequence in read.
- Trim low-quality bases in reads.

## Output:
- `{sample}.clean.fq`: Clean fastq file after trimming.
- `{sample}_summary.txt`: Trim summary.


## Arguments
`--fq` Required. Fastq data.

`--adaptor_3` 3' end adaptor sequence. Default `GCGGAAGCAGTGGTATCAACGCAGAGTACAACAAGGTAC`.

