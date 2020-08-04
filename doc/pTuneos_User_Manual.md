# pTuneos User Manual 


## Table of Contents
1. [General Description](#general-description)  
2. [Dependencies](#dependencies)  
    - [Required software](#required-software)  
    - [Python packages](#python-packages) 
    - [R packages](#r-packages) 
3. [Installation via Docker](#installation-via-docker)  
4. [Installation from source](#installation-from-source)  
5. [Usage](#usage)  
6. [Input Files](#input-files)  
    - [Input Files (WES mode)](#input-files (required))  
    - [Input Files (VCF mode)](#input-files (optional))  
    - [References](#references)
7. [Setting parameters](#seting-parameters) 
8. [Output Files](#output-files)  
    - [Column explanation](#column-explanation)  
9. [Contact](#contact)
10. [Algorithmic Flow Chart](#algorithmic-flow-chart)

## General Description

pTuneos is the state-of-the-art computational pipeline for identifying personalized tumor neoantigens from next-generation sequencing data. With raw whole-exome sequencing data and/or RNA-seq data, pTuneos calculates five important immunogenicity features to construct a machine learning-based classifier (Pre&RecNeo) to predict and prioritize neoantigens recognized by T cell, followed by an efficient score scheme (RefinedNeo) to ealuate naturally processed, MHC presented and T cell recognized probability of a predicted neoepitope.

## Dependencies  

#### Hardware:
pTuneos currently tested on x86_64 on ubuntu 16.04.

Note: if you want to call neoantigen from raw sequence data (WES mode), please ensure that you server have enough CPU core and RAM (We tested on machine with 2 CPU, 88 core and 256 GB RAM). 

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
* [sequenza](https://cran.r-project.org/web/packages/sequenza/index.html)
* [squash](https://CRAN.R-project.org/package=squash)


## Installation via Docker
Docker image of pTuneos is at https://cloud.docker.com/u/bm2lab/repository/docker/bm2lab/ptuneos.

1. Install Docker on your computer and make sure it works.

2. Call docker `pull bm2lab/ptuneos` which will download the Docker image.

3. Run the image in interactive mode with your dataset:
        
		docker run -it -v /your/path/to/dataset/:/home/bioworker/dataset bm2lab/ptuneos /bin/bash

4. Change directory into /home/bioworker/project/pTuneos:

		cd /home/bioworker/project/pTuneos

5. Download reference data:

		bash data_download.sh

6. Edit `config_WES.yaml` or `config_VCF.yaml` and fill the proper path of input files.

7. Run the program with follow commands:

		python pTuneos.py WES -i config_WES.yaml

	or

		python pTuneos.py VCF -i config_VCF.yaml


## Installation from source

1. Install all software, python packages and R packages listed above, and make sure each software and package works in your system. 
2. Install multiprocessing and other packages with the `pip` command:

        pip install -U multiprocessing
        pip install -U pyper
        ...
3. Install R package `copynumber` and its dependence:

        source("http://bioconductor.org/biocLite.R")
        biocLite("copynumber")

   Install R package `squash` and `sequenza`:
   
        install.packages('squash')
        install.packages('sequenza')
 
4. Download or clone the pTuneos repository to your local system:

        git clone https://github.com/bm2-lab/pTuneos.git

5. Reference data includes genome fasta, cDNA, peptide(GRCh38 build) could be downloaded and processed through our script.(you should be aware that the version of VEP library you use should match the references used (peptide and cDNA). E.g. if you install VEP release-89, then you should set the `VEP_release` to `release-89`), then you can run:

        bash data_download.sh
        
    a few reference data would be in the fold `database` and processed by custom script in order to run the pipeline, including:

        [Fasta] 
	
        This fold contains the reference fasta file, its bwa index and some other files result from `huamn.fasta`:
        human.fasta	
        human.fasta.amb	
        human.fasta.ann	
        etc...
	
        [VCF_annotation] 
	
        This fold contains the vcf file and its index files used to run GATK best practice:
        dbsnp_138.hg38.vcf.gz
        1000G_phase1.snps.high_confidence.hg38.vcf.gz
        Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
     
        [Protein] 
	
        This fold contains the reference cDNA and protein sequence of human:
        human.cdna.all.fa
        human.pep.all.fa

6. Among the required software listed above, BWA, GATK 3.8, kallisto, picard, samtools, tabix, trimmomatic-0.36, blast and  VarScan.v2.4.2 were prepared in software directory, other software should be installed by user own due to complexity, please refer to the software links above.

7. Fill in the `config_WES.yaml` file with your local path, make sure you have installed all above software and have downloaded reference data.You should be aware that the version of VEP library you use should match the references used (peptide and cDNA). E.g. in the example above used version/release 89 of GRCh38.


## Usage

pTuneos has two modes: `WES` mode and `VCF` mode.

`WES` mode accepts WES and RNA-seq sequencing data as input, it conduct sequencing quality control, mutation calling, hla typing, expression profiling and neoepitope prediction, filtering, annotation.

`VCF` mode accepts mutation VCF file, expression profile, copy number profile and tumor cellularity as input, it performs neoepitope prediction, filtering, annotation directly on input file.

You can use these two modes by:

        python pTuneos.py WES -i config_WES.yaml
        python pTuneos.py VCF -i config_VCF.yaml

## Input Files

### Input Files (WES mode) 
Pair-end matched tumor-normal whole exome sequencing file should be provided for basic neoepitopes identification, expression profile file or raw RNA sequencing file (pairend or single-end) is optional if you want to get expressed neoepitope. pTuneos accepts pair-end matched tumor-normal whole exome sequencing as input. It could be in `.fastq.gz` or `.fastq` format. 
You should specify the right path to the sequencing file in `config_WES.yaml` like:

    #your path to first tumor fastq file
    tumor_fastq_path_first: ~/ncbi/dbGaP-14145/sra/SRR2770550_1.fastq.gz
    #your path to second tumor fastq file
    tumor_fastq_path_second: ~/ncbi/dbGaP-14145/sra/SRR2770550_2.fastq.gz
    #your path to first normal fastq file
    normal_fastq_path_first: ~/ncbi/dbGaP-14145/sra/SRR2669057_1.fastq.gz
    #your path to second normal fastq file
    normal_fastq_path_second: ~/ncbi/dbGaP-14145/sra/SRR2669057_2.fastq.gz
    #your path to first RNA-seq fastq file
    tumor_rna_fastq_1: ~/ncbi/dbGaP-14145/sra/SRR2673065_1.fastq.gz
    #your path to second RNA-seq fastq file
    tumor_rna_fastq_2: ~/ncbi/dbGaP-14145/sra/SRR2673065_2.fastq.gz
We give the downloading script of WES+RNA-seq testing data in fold `WES_example_data/`. You can use `bash test_data_download.sh` to download these data for testing. (Note: If your RNA-seq data was single-end, just set `tumor_rna_fastq_2` to `None`. In addition, if you know the `fragment length` and `Standard deviation of fragment length` of your single end RNA-seq, replace it with your values, otherwise, just leave them unchanged.)

### Input Files (VCF mode)
Input file for `VCF` mode contains:
* mutaiton file in vcf format from mutect2.
* expression profile in the format same as mentioned in WES mode (recommend obtain from `kallisto`).
* copynumber profile (recommend obtain from `sequenza`).
* tumor cellularity (bewteen 0 and 1) (also recommend obtain from `sequenza`).

We give the example data of these files in fold `VCF_example_data/`.

### References 
The following references are required for pTuneos to run:
* Reference DNA sequence and its annotation file. These files are used in somatic variant calling process.

        [Genome reference]
        human.fasta
        dbsnp_138.hg38.vcf.gz
        1000G_phase1.snps.high_confidence.hg38.vcf.gz
        Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
        CosmicCodingMuts_chr_M_sorted.vcf.gz 
    
* Peptide and cDNA: The peptide reference is a FASTA file containing all peptides and all cDNA sequences of the human proteome.

        [cDNA and protein]
        Homo_sapiens.GRCh38.cdna.all.fa
        Homo_sapiens.GRCh38.pep.all.fa
        
* EnsemblVEP: VEP cache database (It should be emphasized that it is of very high importance that the references and VEP
match in release version (e.g. release-89)).

        [EnsemblVEP]
        homo_sapiens_vep_89_GRCh38.tar.gz


## Setting parameters
User should set all the parameters in the configuration file `config_WES.yaml` or `config_VCF.yaml`. The configuration file contains three parts of parameters:


* Input data parameters, including path of DNA/RNA sequencing data, output fold, run name, hla alleles, expression file and thread number (for WES mode). If your RNA-seq data was single-end, just set `tumor_rna_fastq_2` to `None`. In addition, if you know the `fragment length` and `Standard deviation of fragment length` of your single end RNA-seq, replace it with your values, otherwise, just leave them unchanged.
(Note: user could specific hla allele throught `hla_str`, otherwise set it to `None`, the pipeline will make the prediction utilizing sequencing data. If RNA sequencing data is provided, please also set expression file to `None`.)
* Some filter parameter including mutation sequence depth, mutation variant allele fraction(vaf), binding affinity rank and expression FPKM.
* Software excutable path of opitype, vep, netMHCpan and PyClone.

## Output Files 
pTuneos output four result files contain information of identified neoepitopes corresponding to nonsynonymous point mutation and INDEL mutation.

The output files are the following: 
final_neo_model.tsv 

    The file is a TSV file with the extracted mutated peptides derived from nonsynonymous point mutation and INDEL mutation with a model-based
    score measures the immunity of neoepitopes.



### Column explanation

The prediction output (final_neo_model.tsv) for each peptide pair consists of the following columns:

| Column Name           | Description |
| -----------           | ----------- |
| Position              | Mutation position in genome. |
| HLA_type              | HLA allele name. |
| Gene                  | HUGO symbol name of mutatied gene. |
| WT_pep                | The extracted normal peptide. |
| WT_Binding_EL         | %Rank of prediction score for nomal peptides use NetMHCpan4.0 (defalut model). |
| WT_Binding_Rank       | %Rank of prediction score for nomal peptides use NetMHCpan4.0 (-ba model). |
| MT_pep                | The extracted mutant peptide. |
| MT_Binding_EL         | %Rank of prediction score for mutated peptides use NetMHCpan4.0(defalut model). |
| MT_Binding_Rank       | %Rank of prediction score for mutant peptides use NetMHCpan4.0 (-ba model). |
| Transcript_name       | Ensembl transcript ID |
| Mutation              | Necleotide change of mutated gene |
| AA_change             | Amino acid change annotated in VEP file. |
| Variant_allele_frequency  | Genomic allele frequency detected by MuTect2. |
| DriverGene_Lable      | TRUE if the HUGO symbol is in the cosmic reference list, FALSE if it is not. |
| MT_Binding_level_des  | Binding level description of mutated peptide. |
| WT_Binding_level_des  | Binding level description of normal peptide. |
| Homolog_pep           | The extracted homologous peptide of neo-peptide in human protein. |
| Homolog_Binding_EL    | %Rank of prediction score for homologous peptides use NetMHCpan4.0 (defalut model).|
| Recognition_score     | T cell recognition score calculated based on TCR cross reactivity. |
| Hydrophobicity_score  | Neo-peptide immunity mesurement based on animo acid hydrophobicity. |
| Self_sequence_similarity | Sequence similarity bewteen mutated peptide and normal(homglogous) peptide, We select the bigger one as final score |
| Model_pro       | Model prediction score (probability) for neoepitope recognized by T cell denpend on Recognition_score, Hydrophobicity_score, Self_sequence_similarity, WT_Binding_EL, MT_Binding_EL (calculated by Pre&RecNeo). |
| Immuno_effect_score	| Refined immunogenicty score for neoepitopes (calculated by RefinedNeo). |


## Contact   
 
1410782Chiz@tongji.edu.cn or qiliu@tongji.edu.cn

Biological and Medical Big data Mining Lab  
Tongji University    

## Algorithmic Flow Chart

![](pTuneos_flow_chart.jpg)
