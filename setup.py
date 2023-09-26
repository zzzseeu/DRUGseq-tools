import setuptools
from drug.__init__ import __VERSION__, ASSAY_DICT

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as fp:
    install_requires = fp.read()

entrys = ['drug=drug.drug:main',]
# for assay in ASSAY_DICT:
#     entrys.append(f'multi_{assay}=drug.{assay}.multi_{assay}:main')
entry_dict = {
        'console_scripts': entrys,
}


setuptools.setup(
    name="DRUGseq-tools",
    version=__VERSION__,
    author="zhouxin",
    author_email="xinzhou@xremed.com",
    description="Transcriptome Analysis Pipeline for Drug Discovery",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    include_package_data=True,
    entry_points=entry_dict,
    install_requires=install_requires,
)