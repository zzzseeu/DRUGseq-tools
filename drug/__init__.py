import os
__VERSION__ = "0.0.1"
ASSAY_DICT = {"plastech": "Plastech-seq: Transcriptome sequencing analysis based on barcode and UMI.",
              "phd": "PHD-seq: Target transcriptome sequencing analysis based on barcode and UMI."}


RUN_THREADS = {
    'trim': 5,
    'mapping': 20,
    'count': 5
}

ROOT_DIR = os.path.dirname(__file__)