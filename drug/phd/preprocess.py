import gzip
import statistics
from collections import Counter, defaultdict

import pandas as pd
import pysam

from drug.phd.__init__ import __STEPS__
from drug.toolkits.report import reporter
from drug.toolkits.utils import *

PRIMER_LEN = 17

def correct_seq(error_seq: str, seq_dict: dict) -> str:
    """
    ## Description:
        Replace incorrect barcode sequence with right barcode.
        
    ## Parameters:
        error_seq (str): incorrect barcode sequence with mismatch base pairs.
        
    ## Return:
        raw_seq (str): correct barcode sequence.
    """
    raw_seq = seq_dict[error_seq]
    return raw_seq


class PREPROCESS:
    """
    Features: 
    - Extract the Barcode and UMI information in R1, and use it as the header of R2 read.
    - Filter barcode: Only one base mismatch is allowed at most, and the mismatched base must be a low-quality base.


    Outputs:
    - `{sample}.fq(.gz)` R2 data with modified read header.
    - `stat.txt` Barcode summary.
    """
    def __init__(self, args, step):
        self.step = step
        # required parameters
        self.plate_tsv = args.plate_tsv
        self.sample = args.sample
        # default parameters
        self.outdir = args.outdir
        self.barcode_range = args.barcode_range.split(',')
        self.n_mismatch = 1
        self.min_qual = int(args.min_qual)
        self.gzip = args.gzip
        
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
        df['Tag'] = df['Plate'].astype(str) + '-' + df['Well'].astype(str)
        sample_dict = df.to_dict('list')
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
                          read='1')
            f2 = get_read(library_id=l,
                          library_path=p,
                          read='2')
            if len(f1) != len(f2):
                raise Exception(f"{l} Read1 and Read2 fastq number do not match!")
            fastq_dict[l]['fq1']=f1
            fastq_dict[l]['fq2']=f2
            fastq_dict[l]['info'] = info
        
        return fastq_dict
          
          
    @logit    
    def run(self):
        # prepare
        barcode_range = self.barcode_range
        fastq_dict = self.parse_fastq(plate_tsv=self.plate_tsv)
        # output
        # check outdir
        check_dir(f'{self.out_prefix}')
        # store fastq file (gzip)
        if self.gzip:
            out_fq = gzip.open(f'{self.out_prefix}/{self.sample}.fq.gz', "wt")
        else:
            out_fq = open(f'{self.out_prefix}/{self.sample}.fq', "wt")
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
                    PREPROCESS.run.logger.info(f'Processed {tmp} reads.')
                f1_seq = entry1.sequence
                f1_qual = entry1.quality
                barcode = f1_seq[int(barcode_range[0])-1:int(barcode_range[1])]
                bc_qual = f1_qual[int(barcode_range[0])-1:int(barcode_range[1])]
                
                ### At most one base mismatch is allowed, and the base must be a low-quality base.
                # filter barcode:
                # 1. only barcode in the barcode dict is accepted.
                # 2. only one missmatch is allowed in the barcode.
                # 3. the missmatch base is a low quality base.
                if barcode in barcode_dict:
                    bc_qual = [ord(i)-33 for i in bc_qual]
                    bc_q30 = sum(bc_qual)/len(bc_qual)
                    diff_idx = [i for i in range(len(barcode)) if barcode[i]!=barcode_dict[barcode][i]]

                    if diff_idx==[]:
                        if bc_q30>=30:
                            s['barcode Q30'] += 1
                    elif diff_idx!=[]:
                        if bc_qual[diff_idx[0]]<10:
                            if bc_q30>=30:
                                s['barcode Q30'] += 1                      
                        else:
                            continue
                    s['valid reads'] += 1
                    reads_counter[barcode2tag[barcode_dict[barcode]]] += 1
                    # write valid reads
                    new_head = f'@{barcode2tag[barcode_dict[barcode]]}_{tmp}'
                    new_seq = f'{new_head}\n{entry1.sequence[int(barcode_range[1])+PRIMER_LEN:]}\n+\n{entry1.quality[int(barcode_range[1])+PRIMER_LEN:]}\n'  
                    out_fq.write(f'{new_seq}')      
        
        out_fq.close()
        
        reads_counter = sorted(reads_counter.items(), key = lambda i: -i[1])
                
        ### sum barcode step:
        preprocess_summary = defaultdict()
        preprocess_summary['Total reads:'] = format_int(s['total_reads'])
        valid_reads = s['valid reads']
        valid_reads_percent = round(valid_reads*100/s['total_reads'], 2)
        preprocess_summary['Valid reads:'] = f'{format_int(valid_reads)} ({valid_reads_percent}%)'
        preprocess_summary['Median read counts for barcode:'] = f'{format_int(int(statistics.median([i[1] for i in reads_counter])))}'
        barcode_q30_reads = s['barcode Q30']
        barcode_q30_reads_percent = round(s['barcode Q30']*100/s['total_reads'], 2)
        preprocess_summary['Barcodes Q30:'] = f'{format_int(barcode_q30_reads)} ({barcode_q30_reads_percent}%)'
        preprocess_summary = pd.DataFrame.from_dict(preprocess_summary, orient="index")
        preprocess_summary.to_csv(f'{self.out_prefix}/{self.sample}_summary.txt', sep='\t', header=False)

        preprocess_plot = {'count_values': [i[1] for i in reads_counter],
                        'count_labels': [i[0] for i in reads_counter]}
        
        report = reporter(assay='phd',
                        name=self.step,
                        outdir=f'{self.outdir}',
                        sample=self.sample,
                        stat_file=f'{self.out_prefix}/{self.sample}_summary.txt',
                        plot=preprocess_plot)
        report.get_report()

def preprocess(args):
    step = "preprocess"
    preprocess_obj = PREPROCESS(args, step)
    preprocess_obj.run()
    
    
def get_preprocess_para(parser, optional=False):
    parser.add_argument("--plate_tsv", help="Required. A record file, 3 columns, seperated by tab: \
        Library\tPath\tInfo", required=True)
    parser.add_argument("--barcode_range", help="Barcode position in Read 1. Default `1, 7`.",
                        default="1,7")
    if optional:
        parser.add_argument("--gzip", help="Output gzip fastq file.", 
                            action="store_true")
        parser.add_argument("--min_qual", help="Min quality for barcode base pair. Default `10`.", 
                            default=10)
        parser = common_args(parser)
    return parser