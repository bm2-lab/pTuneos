# pTuneos: prioritizing Tumor neoantigen from next-generation sequencing data #

pTuneos is the state-of-the-art computational pipeline for identifying personalized tumor neoantigens from next-generation sequencing data. With raw whole-exome sequencing data and/or RNA-seq data, pTuneos calculates five important immunogenicity features to construct a machine learning-based classifier (Pre&RecNeo) to predict and prioritize neoantigens recognized by T cell, followed by an efficient score scheme (RefinedNeo) to ealuate naturally processed, MHC presented and T cell recognized probability of a predicted neoepitope.

#### Authors:
Chi Zhou and Qi Liu

#### Citation:
Zhou, C., Wei, Z., Zhang, Z. et al. pTuneos: prioritizing tumor neoantigens from next-generation sequencing data. Genome Med 11, 67 (2019) doi:10.1186/s13073-019-0679-x
#### Web sever:
TBD

## Dependencies

#### Hardware:
pTuneos currently test on x86_64 on ubuntu 16.04.

#### Required software:
* [Python 2.7](https://www.python.org/downloads/release/python-2712/)
* [R 3.2.3](https://cran.r-project.org/src/base/R-3/R-3.2.3.tar.gz)
* [NetMHCpan 4.0](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan)
* [Variant Effect Predictor (VEP)](https://github.com/Ensembl/ensembl-vep)
* [BWA](https://github.com/lh3/bwa)
* [samtools](https://github.com/samtools)
* [Optitype](https://github.com/FRED-2/OptiType)
* [Pyclone](https://bitbucket.org/aroth85/pyclone/wiki/Tutorial)
* [GATK 3.8](https://software.broadinstitute.org/gatk/best-practices/)
* [Picard tools](https://broadinstitute.github.io/picard/)
* [Java 8](https://java.com/en/download/help/linux_x64rpm_install.xml)
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
* [sequenza 2.1.2](https://cran.r-project.org/src/contrib/Archive/sequenza/sequenza_2.1.2.tar.gz)
* [squash](https://CRAN.R-project.org/package=squash)


## Installation

### Install via Docker
Docker image of pTuneos is at https://cloud.docker.com/u/bm2lab/repository/docker/bm2lab/ptuneos.
See the [user manual](/doc/pTuneos_User_Manual.md) for a detailed description usage.

### Install from source
1. Install all software listed above.

2. Download or clone the pTuneos repository to your local system:

        git clone https://github.com/bm2-lab/pTuneos.git

3. Obtain the reference files from GRCh38. These include cDNA, peptide; please refer to [user manual](/doc/pTuneos_User_Manual.md) for a detailed description.


## Usage
pTuneos has two modes, `WES` mode and `VCF` mode.

`PairMatchDna` mode accepts WES and RNA-seq sequencing data as input, it conduct sequencing quality control, mutation calling, hla typing, expression profiling and neoantigen prediction, filtering, annotation.

`VCF` mode accepts mutation VCF file, expression profile, copy number profile and tumor cellularity as input, it performs neoantigen prediction, filtering, annotation directly on input file.

You can use these two mode by:

        python pTuneos.py WES -i config_WES.yaml

or

        python pTuneos.py VCF -i config_VCF.yaml
        
## User Manual 
For detailed information about usage, input and output files, test examples and data
preparation please refer to the [pTuneos User Manual](/doc/pTuneos_User_Manual.md)


## Contact   

1410782Chiz@tongji.edu.cn or qiliu@tongji.edu.cn
Tongji University, Shanghai, China
