# iTunes: identification of personalized Tumor neoantigens from next-generation sequencing data #

iTunes is the state-of-the-art computational pipeline for identifying personalized tumor neoantigens from next-generation sequencing data. With raw whole-exome sequencing data and/or RNA-seq data, iTunes calculates five important immunogenicity features to construct a machine learning-based classifier (vitroneo) to predict and prioritize neoantigens with strong in vitro immunologic effects, followed by an efficient score scheme (vivoneo) to identify neoantigens with in vivo immunologic effects.

#### Authors:
Chi Zhou and Qi Liu

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
* [NetMHCpan 4.0](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan)
* [Variant Effect Predictor (VEP)](https://github.com/Ensembl/ensembl-vep)
* [BWA](https://github.com/lh3/bwa)
* [samtools](https://github.com/samtools)
* [strelka](https://github.com/Illumina/strelka)
* [Optitype](https://github.com/FRED-2/OptiType)
* [Pyclone](https://bitbucket.org/aroth85/pyclone/wiki/Tutorial)
* [GATK 3.8](https://software.broadinstitute.org/gatk/best-practices/)
* [Picard tools](https://broadinstitute.github.io/picard/)
* [Java 8](https://java.com/en/download/help/linux_x64rpm_install.xml)
* [Varscan2](http://varscan.sourceforge.net/)
* [kallisto](http://pachterlab.github.io/kallisto/)
* [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [vcftools](http://vcftools.sourceforge.net/)
* [blast](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* [tabix](http://www.htslib.org/doc/tabix.html)
* [gawk]()

#### Required Python package:
* [yaml](https://pypi.org/project/yaml-1.3/)
* [XGboost](https://pypi.org/project/xgboost/)
* [biopython](https://pypi.org/project/biopython/)
* [scikit-learn==0.19.1](https://pypi.org/project/scikit-learn/)
* [pandas](https://pypi.org/project/pandas/)
* [numpy](https://pypi.org/project/numpy/)
* [imblearn](https://pypi.org/project/imblearn/)
* [Pyomo](https://pypi.org/project/Pyomo/)
* [tables](https://pypi.org/project/tables/)
* [pysam](https://pypi.org/project/pysam/)
* [PypeR](https://pypi.org/project/PypeR/)
* [multiprocessing](https://pypi.org/project/multiprocessing/)
* [subprocess](https://pypi.org/project/subprocess/)
* [math](https://pypi.org/project/math/)
* [matplotlib](https://pypi.org/project/matplotlib/)
* [collections](https://pypi.org/project/collections/)


#### Required R package:
* [cpynumber](http://www.bioconductor.org/packages/release/bioc/html/copynumber.html)
* [sequenza](https://cran.r-project.org/web/packages/sequenza/index.html)
* [squash](https://CRAN.R-project.org/package=squash)


## Installation

1. Install all software listed above.

2. Download or clone the iTunes repository to your local system

        git clone https://github.com/bm2-lab/iTunes.git

3. Obtain the reference files from GRCh38. These include cDNA, peptide; see the References section in the [user manual](/doc/iTunes_User_Manual.md)
for a detailed description.


## Usage
iTunes has two modes, `WES` mode and `VCF` mode.

`PairMatchDna` mode accepts WES and RNA-seq sequencing data as input, it conduct sequencing quality control, mutation calling, hla typing, expression profiling and neoantigen prediction, filtering, annotation.

`VCF` mode accepts mutation VCF file, expression profile, copy number profile and tumor cellularity as input, it performs neoantigen prediction, filtering, annotation directly on input file.

You can use these two mode by:

        python iTunes.py WES -i config_WES.yaml
        python iTunes.py VCF -i config_VCF.yaml
        
## User Manual 
For detailed information about usage, input and output files, test examples and data
preparation please refer to the [iTunes User Manual](/doc/iTunes_User_Manual.md)


## Contact   

1410782Chiz@tongji.edu.cn or qiliu@tongji.edu.cn
Tongji University, Shanghai, China
