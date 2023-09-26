import statistics
import subprocess
from collections import defaultdict

import pandas as pd
import pysam

from drug.phd.__init__ import __STEPS__
from drug.toolkits.report import reporter
from drug.toolkits.utils import *

UMI_LEN = 4

class COUNT:
    """
    Features:
    - Count umi for each gene in each barcode.
    - Filter UMI: 
        1. Cannot contain 'N'.
        2. Cannot be a multimer, such as 'AAAAAAAAAA'.
        3. Cannot have base quality lower than 10.
    """
    def __init__(self, step, args) -> None:
        
        self.step = step
        self.sample = args.sample
        self.outdir = args.outdir
        self.probe = args.probe
        self.fq = args.fq
        self.n_mismatch = int(args.n_mismatch)
        self.clean = args.clean
        
        self.count_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int))) 
        self.unique_mapped_reads = 0 
        self.umi_q30 = 0
        self.valid_mapped_reads = 0
        self.summary_dict = dict()
        
        self.out_prefix = f'{self.outdir}/0{__STEPS__.index(self.step)}.{self.step}'
        check_dir(f'{self.out_prefix}')
        self.fileprefix = f'{self.out_prefix}/{self.sample}'
        
        self.probe_fa = f'{self.fileprefix}.probe.fa'
        self.sam = f'{self.fileprefix}.sam'
        self.bam = f'{self.fileprefix}.bam'
        self.fq_fa = f'{self.fileprefix}.fq.fa'
        
        self.count_detail_file = f'{self.fileprefix}.detail.txt'
        self.count_summary = f'{self.fileprefix}.metadata.txt'
       
    @logit    
    def probe2fa(self) -> None:
        with open(self.probe, 'rt') as f, open(self.probe_fa, 'wt') as ff:
            lines = f.readlines()
            for line in lines:
                s = line.strip('\n').split('\t')
                t = s[1]
                g = s[0]
                ff.write(f'>{g}_{len(t)}\n{t}\n')
                
    @logit         
    def build_index(self) -> None:
        cmd = f'hisat2-build -p 12 -q {self.probe_fa} {self.fileprefix}'
        subprocess.check_call(cmd, shell=True)
        COUNT.build_index.logger.info(cmd)
        
    @logit   
    def mapping(self) -> None:
        cmd = (f'hisat2 -x {self.fileprefix} '
               f'-U {self.fq} '
               f'-S {self.sam} '
               f'-p 15 '
               f'-5 {UMI_LEN} '
               f'-3 68 --quiet --summary-file {self.fileprefix}.align_summary.txt '
               f'--sensitive -k 1 --no-softclip ')
        COUNT.mapping.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
        
    
    @logit
    def sam2bam(self) -> None:
        cmd = (f'samtools view -b {self.sam} > {self.bam}')
        COUNT.sam2bam.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
        
        
    @logit 
    def fastq2fasta(self) -> None:
        fq = pysam.FastxFile(self.fq)
        with open(self.fq_fa, 'wt') as f:
            for entry in fq:
                name = entry.name
                quality = entry.quality
                seq = entry.sequence
                f.write(f'>{name}\n{seq}_{quality}\n')
        
    @logit
    def parse_umi(self) -> None:
        fa = pysam.FastaFile(self.fq_fa)
        bam = pysam.AlignmentFile(self.bam, 'rb')

        for entry in bam:
            name = entry.query_name
            mapping_qual = entry.mapping_quality
            bc = name.split('_')[0]
            attr = fa.fetch(name).strip('\n').split('_')
            seq = attr[0]
            qual = attr[1]
            ref = entry.reference_name
            # if read is not mapped, or the mapping quality lower than 30, or 
            # the mismatch number is more than 4
            if entry.has_tag('NH'):
                if entry.get_tag('NH') > 1:
                    continue
                elif mapping_qual < 60:
                    continue
                elif entry.get_tag('NM') > self.n_mismatch:
                    continue
                else:
                    self.unique_mapped_reads += 1               
                    probe_len = int(ref.split('_')[-1])
                    umi_seq1 = seq[:UMI_LEN]
                    umi_qual1 = qual[:UMI_LEN]
                    umi_seq2 = seq[UMI_LEN+probe_len:UMI_LEN+probe_len+UMI_LEN]
                    umi_qual2 = qual[UMI_LEN+probe_len:UMI_LEN+probe_len+UMI_LEN]
                    umi = umi_seq1+umi_seq2
                    umi_qual = umi_qual1+umi_qual2
                    umi_qual = [ord(i)-33 for i in umi_qual]
                    # filter UMI: 
                    # 1. Must not be a homopolymer, e.g. AAAAAAAAAA
                    # 2. Must not contain N
                    # 3. Must not contain bases with base quality < 10
                    if len(set(umi)) == 1 or 'N' in umi or min(umi_qual)<10:
                        continue
                    if statistics.mean(umi_qual)>=30:
                        self.umi_q30 += 1
                    self.valid_mapped_reads += 1
                    self.count_dict[bc][ref][umi] += 1
        
        self.summary_dict['Unique mapped reads:'] = format_int(self.unique_mapped_reads)
        self.summary_dict['Unique mapped reads with valid UMI:'] = f'{format_int(self.valid_mapped_reads)} ({round(self.valid_mapped_reads*100/self.unique_mapped_reads, 2)}%)'
        self.summary_dict['UMIs Q30:'] = f'{format_int(self.umi_q30)} ({round(self.umi_q30*100/self.unique_mapped_reads, 2)}%)'
                                        
    @logit  
    def count_umi(self) -> None:
        for bc in self.count_dict:
            for g in self.count_dict[bc]:
                correct_umi(self.count_dict[bc][g])
                
        with open(self.count_detail_file, 'wt') as ff:
            ff.write(f'Barcode\tGene_name\tUMI\tCount\n')
            for bc in self.count_dict:
                for g in self.count_dict[bc]:
                    for u in self.count_dict[bc][g]:
                        ff.write(f'{bc}\t{g}\t{u}\t{self.count_dict[bc][g][u]}\n')

                    
    @staticmethod
    def get_df_sum(df, col: str ='UMI') -> None:
        def num_gt2(x):
            return pd.Series.sum(x[x > 1])

        df_sum = df.groupby('Barcode', as_index=False).agg({
            'Count': ['sum', num_gt2],
            'UMI': 'count',
            'Gene_name': 'nunique'
        })
        df_sum.columns = ['Barcode', 'readcount', 'UMI2', 'UMI', 'GeneCount']
        df_sum = df_sum.sort_values(col, ascending=False)
        return df_sum
    
    
    @logit
    def write_matrix(self, df: pd.DataFrame) -> None:
        # output count matrix and count summary
        df_UMI = df.groupby(['Barcode', 'Gene_name'],
                            as_index=False).agg({'UMI': 'count'})
        data = df_UMI.pivot(values='UMI',
                           columns='Barcode',
                           index='Gene_name',).fillna(0).astype(int)
        
        data.to_csv(path_or_buf=f'{self.fileprefix}.matrix.txt',
                sep='\t',
                header =True,
                index=True)
            
    def clean_files(self):
        cmd = (
            f'rm -rf {self.sam} {self.fq_fa} {self.probe_fa} '
            f'{self.fileprefix}.*ht*'
        )
        subprocess.check_call(cmd, shell=True)
        
    @logit
    def run(self):
        self.probe2fa()
        self.build_index()
        self.mapping()
        self.sam2bam()
        self.fastq2fasta()
        self.parse_umi()
        self.count_umi()
        df = pd.read_csv(self.count_detail_file, sep='\t')
        self.write_matrix(df)
        df_sum = self.get_df_sum(df, col='readcount')
        df_sum.to_csv(self.count_summary, sep='\t', index=False)
        
        count_summary = pd.DataFrame.from_dict(self.summary_dict, orient="index")
        count_summary.to_csv(f'{self.out_prefix}/{self.sample}_summary.txt', 
                             sep='\t', header=False)
        report = reporter(assay='phd',
                        name=self.step,
                        outdir=f'{self.outdir}',
                        sample=self.sample,
                        stat_file=f'{self.out_prefix}/{self.sample}_summary.txt')
        report.get_report()
        
        if self.clean:
            self.clean_files()
        
def count(args):
    step = 'count'
    count_obj = COUNT(step, args)
    count_obj.run()
    
    
def get_count_para(parser, optional=False):
    parser.add_argument("--probe", help="Required. Probe sequence file path, 2 columns, gene_name\tprobe_sequence.",
                        required=True)
    parser.add_argument("--fq", help="Required. Fastq file path.",
                        required=True)
    parser.add_argument("--n_mismatch", help="The maximum allowed mismatch of bases. Default `4`.",
                        default=4)
    parser.add_argument("--clean", help="Remove intermediate files.",
                        action="store_true")
    
    if optional:
        parser = common_args(parser)
    return (parser)