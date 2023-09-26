## Features:
- Mapping reads to genome and sort bam file.

## Outputs:
- `{sample}_Aligned.out.bam`: Unsorted bam file.
- `{sample}_Aligned.sortedByCoord.out.bam`: Sorted bam file.
- `{sample}_Log.final.out`: STAR summary file.
- `{sample}_Log.out`: STAR log.
- `{sample}_Log.progress.out`: STAR log.
- `{sample}_SJ.out.tab`:
- `{sample}_summary.txt`: Mapping summary file.


## Arguments
`--clean_fq` Clean R2 fastq file. Required.

`--genomeDir` Required. Genome directory.

