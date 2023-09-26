import glob
import importlib
import logging
import os
import re
import sys
import time
from collections import Counter, defaultdict
from datetime import timedelta
from functools import wraps
from itertools import chain, combinations, product


def correct_umi(umi_dict, percent=0.1):
    """
    Description: Correct umi_dict in place.
    
    Args:
        umi_dict: {umi_seq: umi_count}
        percent: if hamming_distance(low_seq, high_seq) == 1 and
            low_count / high_count < percent, merge low to high.
    Returns:
        n_corrected_umi: int
        n_corrected_read: int
    """
    n_corrected_umi = 0
    n_corrected_read = 0

    # sort by value(UMI count) first, then key(UMI sequence)
    umi_arr = sorted(
        umi_dict.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
    while True:
        # break when only highest in umi_arr
        if len(umi_arr) == 1:
            break
        umi_low = umi_arr.pop()
        low_seq = umi_low[0]
        low_count = umi_low[1]

        for umi_kv in umi_arr:
            high_seq = umi_kv[0]
            high_count = umi_kv[1]
            if float(low_count / high_count) > percent:
                break
            if hamming_distance(low_seq, high_seq) == 1:
                n_low = umi_dict[low_seq]
                n_corrected_umi += 1
                n_corrected_read += n_low
                # merge
                umi_dict[high_seq] += n_low
                del (umi_dict[low_seq])
                break
    return n_corrected_umi, n_corrected_read


def generate_seq_dict(seq_list: list, n_mismatch: int) -> dict:
    """
    ## Description:
        Generate barcode dictionary.
        
    ## Parameters:
        seq_list (list): str list, recording barcode information. 
        n_mismatch (int): maxium allowed mismatch num.
    
    ## Return:
        seq_dict: {mismatch_seq: raw_seq}      
    """

    seq_dict = {}
    for bc in seq_list:
        mismatch_seqs = findall_mismatch(bc, n_mismatch)
        for i in mismatch_seqs:
            seq_dict[i] = bc
    return seq_dict

def get_read(library_id, library_path, read='1'):
    read1_list = [f'_{read}', f'R{read}', f'R{read}_001']
    fq_list = ['fq', 'fastq']
    suffix_list = ["", ".gz"]
    read_pattern_list = [
        f'{library_path}/*{library_id}*{read}.{fq_str}{suffix}' 
        for read in read1_list 
        for fq_str in fq_list 
        for suffix in suffix_list
    ]
    fq_list = [glob.glob(read1_pattern) for read1_pattern in read_pattern_list]
    fq_list = sorted(non_empty for non_empty in fq_list if non_empty)
    fq_list = list(chain(*fq_list))
    if len(fq_list) == 0:
        print("Allowed R1 patterns:")
        for pattern in read_pattern_list:
            print(pattern)
        raise Exception(
            '\n'
            f'Invalid Read{read} path! \n'
            f'library_id: {library_id}\n'
            f'library_path: {library_path}\n'
        )
    return fq_list

def findall_mismatch(seq: str, n_mismatch: int, bases="ACGTN") -> list:
    """
    ## Description:
        choose locations where there's going to be a mismatch using combinations
    and then construct all satisfying lists using product.
    
    ## Parameters:
        seq (str): sequence.
        n_mismatch (int): maxium allowed mismatch num.
        bases (str): default "ACGTN".
        
    ## Return:
        seq_set (list): all mismatch sequence of the raw sequence.
        
    Example:
    >>> seq_set = seq.findall_mismatch("ACG")
    >>> seq_set == answer
    >>> answer = set(["TCG", "AAG", "ACC", "ATG", "ACT", "ACN", "GCG", "ANG", "ACA", "ACG", "CCG", "AGG", "NCG"])
    True
    """
    seq_set = set()
    seq_len = len(seq)
    if n_mismatch > seq_len:
        n_mismatch = seq_len
    for locs in combinations(range(seq_len), n_mismatch):
        seq_locs = [[base] for base in seq]
        for loc in locs:
            seq_locs[loc] = list(bases)
        for poss in product(*seq_locs):
            seq_set.add("".join(poss))
            
    seq_set = list(seq_set)
    return seq_set




def default_args(parser):
    pass

def format_int(n):
    n1 = format(n, ',d')
    return(n1)

def logit(func):
    '''
    logging start and done.
    '''
    logging.basicConfig(level=logging.INFO, 
                        stream=sys.stdout,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', 
                        datefmt='%Y/%m/%d %H:%M:%S')
    module = func.__module__
    name = func.__name__
    logger_name = f"{module}.{name}"
    logger = logging.getLogger(logger_name)

    @wraps(func)
    def wrapper(*args, **kwargs):
        if args and hasattr(args[0], 'debug') and args[0].debug:
            logger.setLevel(10) # debug

        logger.info('start...')
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        used = timedelta(seconds=end - start)
        logger.info('done. time used: %s', used)
        return result

    wrapper.logger = logger
    return wrapper 


@logit
def parse_sample(sample_file: str):
    """
    Args:
        sample_file (string): There are four columns in total, which are sample name, fastq file path, experimental processing, and remark.
        
    Return:
        sample dict (dict)
    """
    
    def split_line(line):
        return line.rstrip('\n').split('\t')
    
    dic = defaultdict(dict)
    
    # check sample file path and format
    try:
        sample_lines = open(sample_file) 
    except IOError:
        print(f'ERROR: No such file: "{sample_file}"!')
    else:
        if not sample_file.endswith('txt'):
            suffix_c = sample_file.split('.')[-1]
            raise Exception(f'ERROR: Invalid file format:.{suffix_c}! .txt file is required for sample file!')
        else:
            line = sample_lines.readline()
            s = split_line(line)
            if len(s) != 3 and len(s) != 4:
                raise Exception(f'ERROR: Invaild separation! Sample file should be separated by tab.')
            while line:
                line = sample_lines.readline()
                s = split_line(line)
                sample = s[0]
                fq1 = s[1].split(',')[0]
                fq2 = s[1].split(',')[1]
                if not os.path.exists(fq1):
                    raise IOError(f'ERROR: No such fastq file: "{fq1}"')
                elif not os.path.exists(fq2):
                    raise IOError(f'ERROR: No such fastq file: "{fq2}"')
                treat = s[2]
                if len(s) == 4:
                    remark = s[3]
                else:
                    remark = None
                dic[sample]['fq1'] = fq1
                dic[sample]['fq2'] = fq2
                dic[sample]['treat'] = treat
                dic[sample]['remark'] = remark
                
                return dic
            

def find_assay_init(assay):
    init_module = importlib.import_module(f"drug.{assay}.__init__")
    return init_module


def find_step_module(assay, step):
    init_module = find_assay_init(assay)
    try:
        step_module = importlib.import_module(f"drug.{assay}.{step}")
    except ModuleNotFoundError:
        try:
            step_module = importlib.import_module(f"drug.toolkits.{step}")
        except ModuleNotFoundError:
            module_path = init_module.IMPORT_DICT[step]
            step_module = importlib.import_module(f"{module_path}.{step}")

    return step_module

         
def common_args(parser):
    ### 这是一些通用的参数。
    parser.add_argument("-s", "--sample", help="Sample name")
    parser.add_argument("-o", "--outdir", help="Output directory.", default='./')
    parser.add_argument("-t", "--thread", help="Number of threads for each step.", default=5)
    return parser


def check_file(file_name: str) -> None:
    """
    Args: 
        file_name: string. Checks if the file exists, 
                  and raise FileNotFoundError if it does not exist.
                  
    Return:
        None
    """
    s = 'No such file(s): '
    if not os.path.exists(file_name):
        raise FileNotFoundError(f'{s} {file_name}')
    else:
        pass

    
def check_dir(dir_name: str) -> None:
    """
    Args: 
        dir_name: string. Checks if the directory exists, 
                  and creates the directory if it does not exist.
                  
    Return:
        None
    """
    if not os.path.exists(dir_name):
        os.system(f'mkdir -p {dir_name}')
            
def iter_readline(file):
    with open(file, 'r') as f:
        while True:
            line = f.readline()
            if line:
                yield line
            else:
                return
            
            
def get_id_name_dict(gtf_file):
    """
    get gene_id:gene_name from gtf file
        - one gene_name with multiple gene_id: "_{count}" will be added to gene_name.
        - one gene_id with multiple gene_name: error.
        - duplicated (gene_name, gene_id): ignore duplicated records and print a warning.
    Returns:
        {gene_id: gene_name} dict
    """

    gene_id_pattern = re.compile(r'gene_id "(\S+)";')
    gene_name_pattern = re.compile(r'gene_name "(\S+)"')
    gene_biotype_pattern = re.compile(r'gene_biotype "(\S+)"')
    id_name = defaultdict(dict)
    c = Counter()
    with open(gtf_file) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('#'):
                continue
            tabs = line.split('\t')
            gtf_type, attributes = tabs[2], tabs[-1]
            if gtf_type == 'gene':
                gene_id = gene_id_pattern.findall(attributes)[-1]
                gene_names = gene_name_pattern.findall(attributes)
                gene_biotype = gene_biotype_pattern.findall(attributes)[-1]
                if not gene_names:
                    gene_name = gene_id 
                else:
                    gene_name = gene_names[-1]
                c[gene_name] += 1
                if c[gene_name] > 1:
                    if gene_id in id_name:
                        assert id_name[gene_id] == gene_name, (
                                'one gene_id with multiple gene_name '
                                f'gene_id: {gene_id}, '
                                f'gene_name this line: {gene_name}'
                                f'gene_name previous line: {id_name[gene_id]}'
                            )
                        get_id_name_dict.logger.warning(
                                'duplicated (gene_id, gene_name)'
                                f'gene_id: {gene_id}, '
                                f'gene_name {gene_name}'
                            )
                        c[gene_name] -= 1
                    else:
                        gene_name = f'{gene_name}_{c[gene_name]}'
                id_name[gene_id]['gene_name'] = gene_name
                id_name[gene_id]['gene_biotype'] = gene_biotype
    return id_name


def hamming_distance(string1, string2):
    distance = 0
    length = len(string1)
    length2 = len(string2)
    if (length != length2):
        print(string1, string2)
        raise Exception(f"string1({length}) and string2({length2}) do not have same length")
    for i in range(length):
        if string1[i] != string2[i]:
            distance += 1
    return distance


def add_dict_item(dic: dict, content: str, value: str) -> dict:
    dic[content] = value
    return dic
        
    
def find_step_module_with_folder(assay, step):
    init_module = find_assay_init(assay)
    folder = ""
    try:
        step_module = importlib.import_module(f"drug.{assay}.{step}")
        folder = assay
    except ModuleNotFoundError:
        try:
            step_module = importlib.import_module(f"drug.toolkits.{step}")
            folder = 'toolkits'
        except ModuleNotFoundError:
            module_path = init_module.IMPORT_DICT[step]
            step_module = importlib.import_module(f"{module_path}.{step}")
            folder = module_path.split('.')[1]

    return step_module, folder