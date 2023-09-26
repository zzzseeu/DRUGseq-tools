# Quick start

## Usage Example

### Plastech-seq

```
micromamba activate xrzs

drug plastech barcode --plate_tsv ./fastq_file.tsv -s sample -o results/

drug plastech trim --fq results/00.barcode/sample.fq -s sample -o results/

drug plastech mapping --clean_fq results/01.trim/sample.clean.fq -s sample -o results/ --genomeDir ~/genome/mouse/

drug plastech featureCounts -s sample -o results/ --gtf ~/genome/mouse/Mus_musculus.GRCm39.104.chr.gtf --input_bam results/02.mapping/sample_Aligned.sortedByCoord.out.bam

drug plastech count --bam results/03.featureCounts/sample_name_sorted.bam --gtf ~/genome/mouse/Mus_musculus.GRCm39.104.chr.gtf -s sample -o results/
```

`--plate_tsv` Required. A record file, 3 columns, seperated by tab: \
        Library\tPath\tInfo.

`--fq` Required. Fastq data.

`--clean_fq` Clean R2 fastq file. Required.

`--genomeDir` Required. Genome directory.

`--gtf` GTF path. Required.

Save scripts above as `{sample}.sh`.

You can start your analysis by running:
```
sh ./{sample}.sh
```

### PHD-seq

```
micromamba activate xrzs

drug phd preprocess --plate_tsv data/fastq_file.tsv -s test -o results/

drug phd count --probe /data/Database/PHD-seq/Probes_sequence_files/006_probes.txt --fq results/00.preprocess/test -s test -o results/
```

`--plate_tsv` Required. A record file, 3 columns, seperated by tab: \
        Library\tPath\tInfo.

`--probe` Required. Probe sequence file path, 2 columns, gene_name\tprobe_sequence.

Save scripts above as `{sample}.sh`.

You can start your analysis by running:
```
sh ./{sample}.sh
```

## Note

Before analysis, you need to prepare two tab-seperated files: plate_tsv and plate_info.

- plate_tsv details: 3 columns, Library\tPath\tInfo.

```
(xrzs) xinzhou@xrzs:~/project/056-pro$ cat fastq_file.tsv
Library Path    Info
056-PT25-P1     /data_DELL_MD1400_storage_202208/xinzhou/project/056-pro/data/056-PT25-P1       /data_DELL_MD1400_storage_202208/xinzhou/project/056-pro/data/056-rj.info.txt

```

- plate_info details: 5 columns, Library\tPlate\tWell\tBarcode\tTreatment.

```
Library Plate   Well    Barcode Treatment
056-PT25-P1     RJ20230718-PTS  A1      AACAAGGTAC      treat
056-PT25-P1     RJ20230718-PTS  A2      AACAATCAGG      treat
056-PT25-P1     RJ20230718-PTS  A3      AACATGGAGA      treat
056-PT25-P1     RJ20230718-PTS  A4      AACATTACCG      treat
```