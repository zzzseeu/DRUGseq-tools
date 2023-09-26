import argparse
import inspect
import os
from collections import defaultdict

import drug.toolkits.utils as utils
from drug.drug import ArgFormatter
from drug.__init__ import ASSAY_DICT

PRE_PROCESSING_STEPS = ('sample')
DOCS_ROOT = 'docs'
MANUAL_MD = f'{DOCS_ROOT}/manual.md'
MANUAL_TEMPLATE = f'{DOCS_ROOT}/manual_template.md'


def generate_single_step_doc(assay, step):
    """
    Returns: 
        - md file relative to DOCS_ROOT 
    """
    step_module, folder = utils.find_step_module_with_folder(assay, step)
    func_opts = getattr(step_module, f"get_{step}_para")
    
    class_docs = get_class_docs(step_module)
    argument_docs = get_argument_docs(func_opts)

    folder_path = f'{DOCS_ROOT}/{folder}/'
    if not os.path.exists(folder_path):
        os.system(f'mkdir -p {folder_path}')

    out_md = f'{DOCS_ROOT}/{folder}/{step}.md'
    with open(out_md, 'w') as out_file:
        out_file.write(class_docs)
        out_file.write(argument_docs)
    return f'{folder}/{step}.md'

def get_argument_docs(func_opts):
    argument_docs = ""
    parser = argparse.ArgumentParser(description='drug',formatter_class=ArgFormatter)
    func_opts(parser, optional=False)
    for argument in parser._option_string_actions:
        if not argument in ['-h', '--help']:
            help_msg = parser._option_string_actions[argument].help
            if help_msg:
                help_msg = help_msg.strip()
            argument_docs += (f'`{argument}` {help_msg}\n\n')
    argument_docs = "## Arguments\n" + argument_docs
    return argument_docs


def get_class_docs(step_module):
    titles = ("Features", "Output")
    class_docs = ""
    for child in inspect.getmembers(step_module, inspect.isclass):
        """Filter out class not defined in step_module"""
        class_obj = child[1]
        if class_obj.__module__ != step_module.__name__:
            continue
        doc = inspect.getdoc(class_obj)
        if doc and "Features" in doc:
            for line in doc.split('\n'):
                for title in titles:
                    if line.find(title) != -1:
                        class_docs += "## " + line + '\n'
                        break
                else:
                    class_docs += line + '\n'
    class_docs += '\n\n'
    return class_docs


def write_step_in_manual(md_path, step, manual_handle):
    """
    - [mkref](rna/mkref.md)
    """
    if not step in PRE_PROCESSING_STEPS:
        manual_handle.write(f'- [{step}]({md_path})\n')
    


@utils.logit
def generate_all_docs():
    md_path_dict = defaultdict(dict)
    for assay in ASSAY_DICT:
        init_module = utils.find_assay_init(assay)
        __STEPS__ = init_module.__STEPS__
        generate_all_docs.logger.info(f"Writing docs {assay} ")
        for step in __STEPS__:
            generate_all_docs.logger.info(f"Writing doc {assay}.{step}")
            md_path = generate_single_step_doc(assay, step)
            md_path_dict[assay][step] = md_path
    return md_path_dict

@utils.logit
def write_manual(md_path_dict):
    with open(MANUAL_MD, 'w') as manual_handle:
        with open(MANUAL_TEMPLATE, 'r') as manual_template:
            manual_handle.write(manual_template.read())
        for assay in ASSAY_DICT:
            init_module = utils.find_assay_init(assay)
            steps = init_module.__STEPS__
            manual_handle.write(f'## {ASSAY_DICT[assay]}\n')
            for step in steps:
                write_manual.logger.info(f"Writing manual {assay}.{step}")
                md_path = md_path_dict[assay][step]
                write_step_in_manual(md_path, step, manual_handle)



if __name__ == "__main__":
    md_path_dict = generate_all_docs()
    write_manual(md_path_dict)