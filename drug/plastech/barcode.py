import gzip
import statistics
from collections import Counter, defaultdict

import pandas as pd
import pysam

from drug.plastech.__init__ import __STEPS__
from drug.toolkits.report import reporter
from drug.toolkits.utils import *


# def correct_seq(error_seq: str, seq_dict: dict) -> str:
#     """
#     ## Description:
#         Replace incorrect barcode sequence with right barcode.
        
#     ## Parameters:
#         error_seq (str): incorrect barcode sequence with mismatch base pairs.
        
#     ## Return:
#         raw_seq (str): correct barcode sequence.
#     """
#     raw_seq = seq_dict[error_seq]
#     return raw_seq


class BARCODE:
    """
    Features: 
    - Extract the Barcode and UMI information in R1, and put them in the header of R2 read.
    - Filter barcode: Only one base mismatch is allowed at most, 
                      and the mismatched base must be a low-quality base.
                      
    Outputs:
    - `{sample}.fq`: R2 fastq data after demultiplexing with barcode and UMI in the seq name. 
        Multi samples will only output one fastq file. 
    - `{sample}_summary.txt`: Barcode step summary file.
    """
    def __init__(self, args, step):
        self.step = step
        # required parameters
        self.plate_tsv = args.plate_tsv
        self.sample = args.sample
        # default parameters
        self.outdir = args.outdir
        self.barcode_range = args.barcode_range.split(',')
        self.umi_range = args.umi_range.split(',')
        self.n_mismatch = 1
        self.min_qual = int(args.min_qual)
        self.gzip = args.gzip
        self.out_fq1 = args.out_fq1
        
        self.out_prefix = f'{self.outdir}/0{__STEPS__.index(self.step)}.{self.step}'
        
        
    @logit
    def parse_plate(self, barcode_tsv: str) -> dict:
        """
        ## Description:
            Parse barcode file and return a dict containning plate barcode info.
            
        ## Parameters:
            barcode_tsv (str): barcode file, 3 columns: Sample\tWell\tBarcode\tSmile.
            
        ## Return:
            barcode_dict (dict): barcode dict. Each barcode corresponds to a well on each plate.
        """
        df = pd.read_csv(filepath_or_buffer=barcode_tsv,
                        sep='\t',
                        header=0)
        df['Tag'] = df['Plate'].astype(dtype=str) + '-' + df['Well'].astype(dtype=str)
        sample_dict = df.to_dict(orient='list')
        return sample_dict
        

    @logit
    def parse_fastq(self, plate_tsv: str) -> dict:
        """
        ## Description:
            Parse sample file and return a dict containning fastq path.
            
        ## Parameters:
            plate_tsv (str): sample file, 2 columns: Sample\tPath\tLibrary.
            
        ## Return:
            fastq_dict (dict): fastq file dict.
        """
        fastq_dict = defaultdict(dict)
        df = pd.read_csv(filepath_or_buffer=plate_tsv,
                        sep='\t',
                        header=0)
        df_dict = df.to_dict('list')
        for i in range(len(df_dict['Library'])):
            l = df_dict['Library'][i]
            p = df_dict['Path'][i]
            info = df_dict['Info'][i]
            f1 = get_read(library_id=l,
                          library_path=p,
                          read=1)
            f2 = get_read(library_id=l,
                          library_path=p,
                          read=2)
            if len(f1) != len(f2):
                raise Exception(f"{l} Read1 and Read2 fastq number do not match!")
            fastq_dict[l]['fq1']=f1
            fastq_dict[l]['fq2']=f2
            fastq_dict[l]['info'] = info
        
        return fastq_dict
          
          
    @logit    
    def run(self):
        # prepare
        barcode_range, umi_range = self.barcode_range, self.umi_range
        fastq_dict = self.parse_fastq(plate_tsv=self.plate_tsv)
        # output
        # check outdir
        check_dir(f'{self.out_prefix}')
        # store fastq file (gzip)
        if self.gzip:
            out_fq2 = gzip.open(f'{self.out_prefix}/{self.sample}_2.fq.gz', "wt")
            if self.out_fq1:
                out_fq1 = gzip.open(f'{self.out_prefix}/{self.sample}_1.fq.gz', "wt")
        else:
            out_fq2 = open(f'{self.out_prefix}/{self.sample}_2.fq', "wt")
            if self.out_fq1:
                out_fq1 = open(f'{self.out_prefix}/{self.sample}_1.fq', "wt")

        # statement 
        # total reads
        # valid reads
        # mean reads per barcode
        # Q30 of Barcodes
        # Q30 of UMIs
        s = Counter()
        reads_counter = defaultdict(int)
        # process
        for f in fastq_dict:
            info = fastq_dict[f]['info']
            fq1 = fastq_dict[f]['fq1'][0]
            fq2 = fastq_dict[f]['fq2'][0]
            print(fq1, fq2)
            # parse barcode file
            plate_info = self.parse_plate(info)
            # generate barcode dict
            barcode_dict = generate_seq_dict(plate_info['Barcode'], self.n_mismatch)
            barcode2tag = {plate_info['Barcode'][i]: plate_info['Tag'][i] for i in range(len(plate_info['Barcode']))}      
            # use pysam to read fastq file
            f1, f2 = pysam.FastxFile(fq1), pysam.FastxFile(fq2)
            # performing reads
            for entry1, entry2 in zip(f1, f2):
                s['total_reads'] += 1
                tmp = s['total_reads']
                if tmp % 5000000 == 0:
                    BARCODE.run.logger.info(f'Processed {tmp} reads.')
                f1_seq = entry1.sequence
                f1_qual = entry1.quality
                barcode = f1_seq[int(barcode_range[0])-1:int(barcode_range[1])]
                bc_qual = f1_qual[int(barcode_range[0])-1:int(barcode_range[1])]
                umi = f1_seq[int(umi_range[0])-1:int(umi_range[1])]
                umi_qual = f1_qual[int(umi_range[0])-1:int(umi_range[1])]
                
                ### At most one base mismatch is allowed, and the base must be a low-quality base.
                # filter barcode:
                # 1. only barcode in the barcode dict is accepted.
                # 2. only one missmatch is allowed in the barcode.
                # 3. the missmatch base is a low quality base.
                if barcode in barcode_dict:
                    bc_qual = [ord(i)-33 for i in bc_qual]
                    bc_q30 = sum(bc_qual)/len(bc_qual)
                    diff_idx = [i for i in range(len(barcode)) if barcode[i]!=barcode_dict[barcode][i]]
                    umi_qual = [ord(i)-33 for i in umi_qual]
                    umi_q30 = sum(umi_qual)/len(umi_qual)
                    umi_pass = True
                    # filter UMI:
                    # Must not be a homopolymer, e.g. AAAAAAAAAA
                    # Must not contain N
                    # Must not contain bases with base quality < 10
                    if min(umi_qual)<10 or ('N' in umi) or len(set(list(umi)))==1:
                        umi_pass = False
                    if diff_idx==[] and umi_pass==True:
                        if bc_q30>=30:
                            s['barcode Q30'] += 1
                        if umi_q30>=30:
                            s['UMI Q30'] += 1
                    elif diff_idx!=[] and umi_pass==True:
                        if bc_qual[diff_idx[0]]<10:
                            if bc_q30>=30:
                                s['barcode Q30'] += 1   
                            if umi_q30>=30:
                                s['UMI Q30'] += 1                     
                        else:
                            continue
                    elif diff_idx!=[] and umi_pass==False:
                        continue
                    s['valid reads'] += 1
                    reads_counter[barcode2tag[barcode_dict[barcode]]] += 1
                    # write valid reads
                    new_head = f'@{barcode2tag[barcode_dict[barcode]]}_{umi}_{tmp}'
                    new_seq = f'{new_head}\n{entry2.sequence}\n+\n{entry2.quality}\n'  
                    out_fq2.write(f'{new_seq}')
                    if self.out_fq1:
                        out_fq1.write(f'{new_head}\n{f1_seq[int(umi_range[1]):]}\n+\n{f1_qual[int(umi_range[1]):]}\n')
        
        out_fq2.close()
        
        reads_counter = sorted(reads_counter.items(), key = lambda i: -i[1])
                
        ### sum barcode step:
        barcode_summary = defaultdict()
        barcode_summary['Total reads:'] = format_int(s['total_reads'])
        valid_reads = s['valid reads']
        valid_reads_percent = round(valid_reads*100/s['total_reads'], 2)
        barcode_summary['Valid reads:'] = f'{format_int(valid_reads)} ({valid_reads_percent}%)'
        barcode_summary['Median read counts for barcode:'] = format_int(int(statistics.median([i[1] for i in reads_counter])))
        barcode_q30_reads = s['barcode Q30']
        barcode_q30_reads_percent = round(s['barcode Q30']*100/s['total_reads'], 2)
        barcode_summary['Barcodes Q30:'] = f'{format_int(barcode_q30_reads)} ({barcode_q30_reads_percent}%)'
        umi_q30_reads = s['UMI Q30']
        umi_q30_reads_percent = round(s['UMI Q30']*100/s['total_reads'], 2)
        barcode_summary['UMIs Q30:'] = f'{format_int(umi_q30_reads)} ({umi_q30_reads_percent}%)'
        barcode_summary = pd.DataFrame.from_dict(barcode_summary, orient="index")
        barcode_summary.to_csv(f'{self.out_prefix}/{self.sample}_summary.txt', sep='\t', header=False)

        barcode_plot = {'count_values': [i[1] for i in reads_counter],
                        'count_labels': [i[0] for i in reads_counter]}
        
        report = reporter(assay='plastech',
                        name=self.step,
                        outdir=f'{self.outdir}',
                        sample=self.sample,
                        stat_file=f'{self.out_prefix}/{self.sample}_summary.txt',
                        plot=barcode_plot)
        report.get_report()

def barcode(args):
    step = "barcode"
    barcode_obj = BARCODE(args, step)
    barcode_obj.run()
    
    
def get_barcode_para(parser, optional=False):
    parser.add_argument("--plate_tsv", help="Required. A record file, 3 columns, seperated by tab: \
        Library\tPath\tInfo.", required=True)
    parser.add_argument("--barcode_range", help="Barcode position in Read 1. Default `1, 10`.",
                        default="1,10")
    parser.add_argument("--umi_range", help="UMI position in Read 1. Default `11, 20`.", 
                        default="11,20")
    if optional:
        parser.add_argument("--gzip", help="Output gzip fastq file.", 
                            action="store_true")
        parser.add_argument("--min_qual", help="Min quality for barcode base pair. Default `10`.", 
                            default=10)
        parser.add_argument("--out_fq1", help="Output R1 fastq file without barcode and UMI.", 
                            action="store_true")
        parser = common_args(parser)
    
    return parser
