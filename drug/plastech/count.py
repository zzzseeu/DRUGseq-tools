from collections import defaultdict
from itertools import groupby

import anndata
import hdf5plugin
import pandas as pd
import pysam
from scipy import sparse

from drug.plastech.__init__ import __STEPS__
from drug.toolkits import utils


class COUNT():
    """
    Features:
    - Count umi for each gene in each barcode.
    - Filter UMI: 
        1. Cannot contain 'N'.
        2. Cannot be a multimer, such as 'AAAAAAAAAA'.
        3. Cannot have base quality lower than 10.

    Outputs:
    - `{sample}.detail.txt`: UMI, read count raw file. 4 columns, Barcode\tgeneID\tUMI\tcount.
    - `{sample}.features.txt`: Gene features record file. 3 columns, gene_name\tgene_id\tgene_biotype.
    - `{samples}.metadata.txt`: Sample meta data. 5 columns, Barcode\treadcount\tUMI2(UMI counts>2)\tUMI\tGeneCount.
    - `{sample}.matrix.txt`: Gene expression matrix. Columns are genes, rows are samples.
    - `{sample}.matrix.h5ad`: Gene expression matrix h5 file for downstream analysis.
    """
    def __init__(self, step, args):

        # init
        self.step = step
        self.sample = args.sample
        self.outdir = args.outdir

        # required parameters
        self.bam = args.bam
        self.gtf = args.gtf
        
        # default parameters
        self.outTable = args.outTable

        # output files
        self.outdir_prefix = f'{self.outdir}/0{__STEPS__.index(self.step)}.{self.step}'
        utils.check_dir(f'{self.outdir_prefix}')
        self.fileprefix = f'{self.outdir_prefix}/{self.sample}'
        self.count_detail_file = f'{self.fileprefix}.detail.txt'
        self.count_summary = f'{self.fileprefix}.metadata.txt'
        self.h5file = f'{self.fileprefix}.matrix.h5ad'
    

    @utils.logit
    def bam2table(self):
        """
        bam to detail table
        must be used on name_sorted bam
        """
        samfile = pysam.AlignmentFile(self.bam, "rb")
        with open(self.count_detail_file, 'wt') as fh1:
            fh1.write('\t'.join(['Barcode', 'geneID', 'UMI', 'count']) + '\n')

            def keyfunc(x):
                return x.query_name.split('_', 1)[0]
            for _, g in groupby(samfile, keyfunc):
                gene_umi_dict = defaultdict(lambda: defaultdict(int))
                for seg in g:
                    (barcode, umi) = seg.query_name.split('_')[:2]
                    if not seg.has_tag('XT'):
                        continue
                    gene_id = seg.get_tag('XT')
                    gene_name = seg.get_tag('GN')
                    gene_biotype = seg.get_tag('GT')
                    gene_umi_dict['&'.join([gene_id, gene_name, gene_biotype])][umi] += 1
                for gene_id in gene_umi_dict:
                    utils.correct_umi(gene_umi_dict[gene_id])

                # output
                for gene_id in gene_umi_dict:
                    for umi in gene_umi_dict[gene_id]:
                        fh1.write('%s\t%s\t%s\t%s\n' % (barcode, gene_id, umi,
                                                        gene_umi_dict[gene_id][umi]))
        samfile.close()
        

    @staticmethod
    def get_df_sum(df, col='UMI'):
        def num_gt2(x):
            return pd.Series.sum(x[x > 1])

        df_sum = df.groupby('Barcode', as_index=False).agg({
            'count': ['sum', num_gt2],
            'UMI': 'count',
            'geneID': 'nunique'
        })
        df_sum.columns = ['Barcode', 'readcount', 'UMI2', 'UMI', 'GeneCount']
        df_sum = df_sum.sort_values(col, ascending=False)
        return df_sum


    @utils.logit
    def write_matrix(self, df):
        # output count matrix and count summary
        df_UMI = df.groupby(['Barcode', 'geneID'],
                            as_index=False).agg({'UMI': 'count'})
        data = df_UMI.pivot(values='UMI',
                            columns='Barcode',
                            index='geneID',).fillna(0).astype(int)
        data.insert(0, 'gene_id', list(map(lambda x: x.strip().split('&')[0], data.index.tolist())))
        data.insert(0, 'gene_name', list(map(lambda x: x.strip().split('&')[1], data.index.tolist())))
        data.insert(0, 'gene_biotype', list(map(lambda x: x.strip().split('&')[2], data.index.tolist())))
        
        data = data.set_index(['gene_name', 'gene_id', 'gene_biotype'])
        mtx = sparse.csr_matrix(data.values.T)
        obs = pd.DataFrame()
        obs.index = data.columns.tolist()
        var = data.index.to_frame(index=False).set_index('gene_name')
        adata = anndata.AnnData(mtx, obs=obs, var=var, dtype='int')
        adata.write_h5ad(filename=self.h5file,
                         compression=hdf5plugin.FILTERS["zstd"])
        if self.outTable:
            features = data.index.to_frame()
            features.to_csv(path_or_buf=f'{self.fileprefix}.features.txt',
                                                    sep='\t',
                                                    header=True,
                                                    index=False)
            data.index = data.index.levels[0]
            data.T.to_csv(path_or_buf=f'{self.fileprefix}.matrix.txt',
                    sep='\t',
                    header =True,
                    index=True)

    @utils.logit
    def run(self):
        self.bam2table()
        df = pd.read_csv(self.count_detail_file, sep='\t')
        self.write_matrix(df)
        df_sum = self.get_df_sum(df, col='readcount')
        df_sum.to_csv(self.count_summary, sep='\t', index=False)


def count(args):
    step = 'count'
    count_obj = COUNT(step, args)
    count_obj.run()


def get_count_para(parser, optional=False):
    parser.add_argument("--bam", help="Sorted featureCounts output bamfile.",
                        required=True)
    parser.add_argument("--gtf", help="GTF file path.",
                        required=True)
    parser.add_argument("--outTable", help="Output gene expression count table.",
                        action="store_true")
    if optional:
        parser = utils.common_args(parser)
    return (parser)
