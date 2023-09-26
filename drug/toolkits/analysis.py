import scanpy as sc
import pandas as pd
import omicverse as ov
import hdf5plugin
from scipy import sparse
import numpy as np
import os
from statsmodels.stats.multitest import multipletests
from drug.toolkits import utils


def parse_format(file: str) -> str:
    suffix = file.strip().split('.')[-1]
    return(suffix)          

def read_csv(file: str) -> sc.AnnData:
    """
    Description:
        Read in csv (comma seperated) file and return a AnnData object \
        for downsream analysis.
    
    Parameters:
        file(str): File path.
        
    Return:
        adata: a AnnData object.
        
    Raise:
        IOError: Raise error when invalid file format.
        FileNotFoundError: File does not exist.
    """
    suffix = parse_format(file=file)
    if not suffix=='csv':
        raise IOError('Invaild file format!')
    elif not os.path.exists(path=file):
        raise FileNotFoundError
    else:
        data = pd.read_csv(filepath_or_buffer=file,
                           sep=',',
                           header=0,
                           index_col=0).T
        mtx = data.values
        obs = pd.DataFrame()
        obs.index = data.index.tolist()
        var = pd.DataFrame()
        var.index = data.columns.tolist()
        adata = sc.AnnData(X=mtx, obs=obs, var=var, dtype='int')
        return adata      
    
def read_table(file: str) -> sc.AnnData:
    """
    Description:
        Read in txt/tsv (table seperated) file and return a AnnData object \
        for downsream analysis.
    
    Parameters:
        file(str): File path.
        
    Return:
        adata: a AnnData object.
        
    Raise:
        IOError: Raise error when invalid file format.
        FileNotFoundError: File does not exist.
    """
    # check file format
    suffix = parse_format(file=file)
    if suffix!='txt' and suffix!='tsv':
        raise IOError('Invaild file format!')
    elif not os.path.exists(path=file):
        raise FileNotFoundError
    else:
        data = pd.read_csv(filepath_or_buffer=file,
                           sep='\t',
                           header=0,
                           index_col=0).T
        mtx = data.values
        obs = pd.DataFrame()
        obs.index = data.index.tolist()
        var = pd.DataFrame()
        var.index = data.columns.tolist()
        adata = sc.AnnData(X=mtx, obs=obs, var=var, dtype='int')
        return adata       

def read_h5ad(file: str) -> sc.AnnData:
    """
    Description:
        Read in h5ad file and return a AnnData object \
        for downsream analysis.
    
    Parameters:
        file(str): File path.
        
    Return:
        adata: a AnnData object.
        
    Raise:
        IOError: Raise error when invalid file format.
        FileNotFoundError: File does not exist.
    """
    # check file format
    suffix = parse_format(file=file)
    if suffix!='h5ad':
        raise IOError('Invaild file format!')
    elif not os.path.exists(path=file):
        raise FileNotFoundError
    else:
        adata = sc.read_h5ad(filename=file)
        return adata
    
def read_coldata(coldata: str) -> pd.DataFrame:
    """
    Description:
        Read in coldata file and return a data frame \
        for downsream analysis.
    
    Parameters:
        coldata(str): File path.
        
    Return:
        coldata(pd.DataFrame): Meta.data contains sample group info.
        
    Raise:
        IOError: Raise error when invalid file format.
        FileNotFoundError: File does not exist.
    """
    # check file format
    suffix = parse_format(coldata)
    if not suffix in ['csv', 'txt', 'tsv']:
        raise IOError('Invaild file format!')
    elif not os.path.exists(coldata):
        raise FileNotFoundError
    else:
        if suffix=='csv':
            coldata = pd.read_csv(filepath_or_buffer=coldata,
                                  sep=',',
                                  header=0)
        elif suffix=='txt' or suffix=='tsv':
            coldata = pd.read_csv(filepath_or_buffer=coldata,
                                    sep='\t',
                                    header=0)
        return coldata

def ttest(counts: pd.DataFrame, 
          coldata: pd.DataFrame, 
          treat: str,
          control: str, 
          lfc_threshold: float=1, 
          p_threshold: float=0.05,
          p_type: str='padj', 
          multipletests_method: str='fdr_bh'):
    from scipy.stats import ttest_ind
    
    data=counts
    group1 = coldata[coldata.treatment==treat].barcode.tolist()
    group2 = coldata[coldata.treatment==control].barcode.tolist()
    g1_mean=data[group1].mean(axis=1)
    g2_mean=data[group2].mean(axis=1)
    g=(g2_mean+g1_mean)/2
    g=g.loc[g>0].min()
    fold=(g1_mean+g)/(g2_mean+g)
    #log2fold=np.log2(fold)
    ttest = ttest_ind(data[group1].T.values, data[group2].T.values)
    pvalue=ttest[1]
    padj = multipletests(np.nan_to_num(np.array(pvalue),0), alpha=0.5, 
                        method=multipletests_method, is_sorted=False, returnsorted=False)
    genearray = np.asarray(pvalue)
    result = pd.DataFrame({'pvalue':genearray,'padj':padj[1],'FoldChange':fold})
    result=result.loc[~result['pvalue'].isnull()]
    result['-log(pvalue)'] = -np.log10(result['pvalue'])
    result['-log(padj)'] = -np.log10(result['padj'])
    result['BaseMean']=(g1_mean+g2_mean)/2
    result['log2(BaseMean)']=np.log2((g1_mean+g2_mean)/2)
    result['log2FC'] = np.log2(result['FoldChange'])
    result['abs(log2FC)'] = abs(np.log2(result['FoldChange']))
    result['size']  =np.abs(result['FoldChange'])/10
    result['sig']='NoDiff'
    result.loc[(result[p_type]<p_threshold) & (result['log2FC']>lfc_threshold),
               'sig']='Up'
    result.loc[(result[p_type]<p_threshold) & (result['log2FC']<-lfc_threshold),
            'sig']='Down'
    
    return(result)

def wilcox(counts: pd.DataFrame, 
          coldata: pd.DataFrame, 
          treat: str,
          control: str, 
          lfc_threshold: float=1, 
          p_threshold: float=0.05,
          p_type: str='padj', 
          multipletests_method: str='fdr_bh'):
    from scipy.stats import ranksums
    
    data=counts
    group1 = coldata[coldata.treatment==treat].barcode.tolist()
    group2 = coldata[coldata.treatment==control].barcode.tolist()
    g1_mean=data[group1].mean(axis=1)
    g2_mean=data[group2].mean(axis=1)
    fold=(g1_mean+0.00001)/(g2_mean+0.00001)
    #log2fold=np.log2(fold)
    wilcox = ranksums(data[group1].T.values, data[group2].T.values)
    pvalue=wilcox[1]
    padj = multipletests(np.nan_to_num(np.array(pvalue),0), alpha=0.5, 
                        method=multipletests_method, is_sorted=False, returnsorted=False)
    genearray = np.asarray(pvalue)
    result = pd.DataFrame({'pvalue':genearray,'padj':padj[1],'FoldChange':fold})
    result=result.loc[~result['pvalue'].isnull()]
    result['-log(pvalue)'] = -np.log10(result['pvalue'])
    result['-log(padj)'] = -np.log10(result['padj'])
    result['BaseMean']=(g1_mean+g2_mean)/2
    result['log2(BaseMean)']=np.log2((g1_mean+g2_mean)/2)
    result['log2FC'] = np.log2(result['FoldChange'])
    result['abs(log2FC)'] = abs(np.log2(result['FoldChange']))
    result['size']  =np.abs(result['FoldChange'])/10
    #result=result[result['padj']<alpha]
    result['sig']='NoDiff'
    result.loc[(result[p_type]<p_threshold) & (result['log2FC']>lfc_threshold),
               'sig']='Up'
    result.loc[(result[p_type]<p_threshold) & (result['log2FC']<-lfc_threshold),
            'sig']='Down'

    return result
    
def deseq2(counts: pd.DataFrame,
           coldata: pd.DataFrame,
           treat: str,
           control: str,
           lfc_threshold: float,
           p_threshold: float,
           p_type='padj'
            ):
    """
    Description:
        Use adata for differential expressed genes analysis and output DEGs.
    
    Parameters:
    
    Return:
    
    Raise
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    dds = DeseqDataSet(counts=counts.T,
                       metadata=coldata,
                       design_factors='treatment',
                       n_cpus=8,)
    dds.deseq2()
    stat_res = DeseqStats(dds, 
                          n_cpus=8, 
                          contrast=['treatment', treat, control])
    stat_res.summary()
    result = stat_res.results_df
    result = result.rename(columns={'log2FoldChange':'log2FC'})
    result['sig']='NoDiff'
    result.loc[(result[p_type]<p_threshold) & (result['log2FC']>lfc_threshold),
               'sig']='Up'
    result.loc[(result[p_type]<p_threshold) & (result['log2FC']<-lfc_threshold),
            'sig']='Down'   
    return result
    
def gsea(genelist,
         gene_sets,
         organism):
    import gseapy as gp
    enr = gp.enrichr(gene_list=genelist, # or "./tests/data/gene_list.txt",
                    gene_sets=gene_sets,
                    organism=organism, # don't forget to set organism to the one you desired! e.g. Yeast
                    )
    result = enr.results
    
    return result

def pywgcna():
    pass

class ANALYSIS:
    def __init__(self, args, step) -> None:
        self.step = step
        self.mtx = args.mtx
        self.metadata = args.metadata
        self.treat = args.treat
        self.control = args.control
        self.lfc_threshold = args.lfc_threshold
        self.pCutoff = args.pCutoff
        self.annotation = args.annotation
        self.outdir = args.outdir
        self.species = args.species
        self.method = args.method
        
        utils.check_dir(self.outdir)
        utils.check_dir(f'{self.outdir}/wgcna')
        
    @utils.logit  
    def run(self) -> None:
        # load gene expression matrix
        data = pd.read_csv(filepath_or_buffer=self.mtx,
                   sep='\t',
                   header=0,
                   index_col=0)
        data = data.drop(['gene_id'], axis=1)
        dds=ov.bulk.pyDEG(data)
        dds.drop_duplicates_index()
        print('drop_duplicates_index success')
        # load sample metadata
        info = pd.read_csv(filepath_or_buffer=self.metadata,
                        sep='\t',
                        index_col=False,
                        header=0)
        # perform deseq2 analysis
        treat_samples = info[info.treatment==self.treat].barcode.tolist()
        control_samples = info[info.treatment==self.control].barcode.tolist()
        if self.method=='DEseq2':
            result=dds.deg_analysis(treat_samples,control_samples,method=self.method)
            dds.normalize()
        else:
            dds.normalize()
            result=dds.deg_analysis(treat_samples,control_samples,method=self.method)
        result = result[~result.padj.isna()]
        result.to_csv(path_or_buf=f'{self.outdir}/diff.genes.csv',
                      sep=',',
                      index=True)
        # set cutoff
        # -1 means automatically calculates
        dds.foldchange_set(fc_threshold=self.lfc_threshold,
                   pval_threshold=self.pCutoff,
                   logp_max=100)
        dds.plot_volcano(title='DEG Analysis',
                 figsize=(10,10),
                 plot_genes_num=8,
                 plot_genes_fontsize=12,
                 legend_bbox=[1.3, 0.5])
        # enrichment analysis
        pathway_dict=ov.utils.geneset_prepare(geneset_path=self.annotation,
                                              organism=self.species)
        rnk=dds.ranking2gsea()
        gsea_obj=ov.bulk.pyGSEA(rnk,pathway_dict)
        enrich_res=gsea_obj.enrichment()
        enrich_res.to_csv(f'{self.outdir}/diff.genes.enrichment.csv', 
                          sep=',',
                          index=True)
        gsea_obj.plot_enrichment(num=10,
                        node_size=[10,20,30],
                        cax_loc=1.2,
                        cax_fontsize=10,
                        fig_title='Enrichment',
                        fig_xlabel='Fractions of genes',
                        figsize=(10,10),
                        cmap='YlGnBu',
                        text_knock=2,
                        text_maxsize=50)
        # wgcna
        wgcna_data = dds.data
        wgcna_data = wgcna_data[wgcna_data.index.isin(result.index)]
        gene_wgcna=ov.bulk.pyWGCNA(wgcna_data,save_path=f'{self.outdir}/wgcna')
        gene_wgcna.calculate_correlation_direct(method='pearson',save=False)
        gene_wgcna.calculate_correlation_indirect(save=False)
        gene_wgcna.calculate_soft_threshold(save=False)
        gene_wgcna.calculate_corr_matrix()
        gene_wgcna.calculate_distance()
        gene_wgcna.calculate_geneTree()
        gene_wgcna.calculate_dynamicMods()
        module=gene_wgcna.calculate_gene_module()
        gene_wgcna.plot_matrix()
        module.to_csv(f'{self.outdir}/module.csv',
                      sep=',',
                      index=True)
        
def analysis(args):
    step = "analysis"
    analysis_obj = ANALYSIS(args, step)
    analysis_obj.run()
    
def get_analysis_para(parser, optional=False):
    parser.add_argument("--mtx", help="Gene expression matrix", 
                        required=True)
    parser.add_argument("--metadata", help="Sample info matrix.",
                        required=True)
    parser.add_argument("--treat", help="Treat group name.", 
                        default="PC")
    parser.add_argument("--control", help="Control group name.", 
                        default="NC")
    parser.add_argument("--annotation", help="Enrichr annotation file", 
                    required=True)
    if optional:
        parser.add_argument("--lfc_threshold", help="Log2FoldChange threshold for DEGs.", 
                            default=1)
        parser.add_argument("--pCutoff", help="pvalue/qvalue threshold for DEGs", 
                            default=0.05)
        parser.add_argument("--species", help="pvalue/qvalue threshold for DEGs", 
                            default='Human')
        parser.add_argument("--method", help="pvalue/qvalue threshold for DEGs", 
                            default='DEseq2')
        parser = utils.common_args(parser)
    
    return parser