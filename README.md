# iTunes: identification of Tumor neoantigens pipeline from next-generation sequencing data #

Given matched tumor-normal whole exome sequencing and tumor RNA-seq sequencing data 
as input, ITUNES infers HLA sub-types, mutated peptides (neo-peptide), variant allele
frequency, expression profile etc feature information. Based on these feature, a model-based
refined ranking-score scheme could identify which of the neo-peptides have strong 
immunogecity. Detailed information please refer to citation.

#### Authors and Email:
Chi Zhou 

#### Citation:
iTunes: identification of personalized Tumor neoantigens from next-generation sequencing data, Submitted, 2018.

#### Web sever:
TBD

## Dependencies

#### Hardware:
iTunes currently test on x86_64 on ubuntu 16.04.

#### Required software:
* [Python 2.7](https://www.python.org/downloads/release/python-2712/)
* [R 3.2.3](https://cran.r-project.org/src/base/R-3/R-3.2.3.tar.gz)
* [NetMHCpan 3.0](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan)
* [Variant Effect Predictor (VEP)](https://github.com/Ensembl/ensembl-vep)
* [bwa](https://github.com/lh3/bwa)
* [samtools](https://github.com/samtools)
* [STAR](https://github.com/alexdobin/STAR)
* [strelka](https://github.com/Illumina/strelka)
* [opitype](https://github.com/FRED-2/OptiType)
* [pyclone](https://bitbucket.org/aroth85/pyclone/wiki/Tutorial)
* [GATK 3.7](https://software.broadinstitute.org/gatk/best-practices/)
* [Picard tools](https://broadinstitute.github.io/picard/)
* [Java 8](https://java.com/en/download/help/linux_x64rpm_install.xml)

#### Required software:
* [yaml]()
* [XGboost]()
* [Bio]()
* [sklearn]()
* [pandas]()
* [numpy]()

## Installation

1. Install all software listed above.

2. Download or clone the iTunes repository to your local system

        git clone https://github.com/bm2-lab/iTunes.git

3. Obtain the reference files from GRCh38. These include cDNA, peptide and COSMIC
files; see the References section in the [user manual](/doc/iTunes_User_Manual.md)
for a detailed description.

4. Fill in the config.yaml file  
    see the config.yaml file for more information

## Usage

Here is a simple example in which somatic mutation calls and gene expression data are
provided, and MHC binding is predicted for HLA types HLA-A01:01 and HLA-B08:01. 

    python path/to/iTunes.py -i config.yaml

All options can be displayed using the usage information with the `-h` option:   

    python path/to/iTunes.py -h


## User Manual 
For detailed information about usage, input and output files, test examples and data
preparation read the [iTunes User Manual](/doc/iTunes_User_Manual.md)


## Contact   

1410782Chiz@tongji.edu.cn or qiliu@tongji.edu.cn
