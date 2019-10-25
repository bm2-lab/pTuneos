#######This script download the testing data for running pTuneos from WES+RNA-seq sequence data###########
#The testing data was from a tumor patient (NCI-3784) of melanoma, which was public available, more detail please refer to 
#Gros et al. Prospective identification of neoantigen-specific lymphocytes in the peripheral blood of melanoma patients. Nature medicine, 2016.

######Please ensure that you have installed sratoolkit, otherwise refer to https://hpc.nih.gov/apps/sratoolkit.html
## tumor exome
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX132/SRX1322495/SRR2602447/SRR2602447.sra
## matched normal exome
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX132/SRX1322495/SRR2602449/SRR2602449.sra
## RNA-seq
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX132/SRX1323066/SRR2603944/SRR2603944.sra

######Extract SRA file into fastq file using fastq-dump, then compress into fastq.gz
fastq-dump -I --origfmt --split-files --gzip SRR2602447.sra
fastq-dump -I --origfmt --split-files --gzip SRR2602449.sra
fastq-dump -I --origfmt --split-files --gzip SRR2603944.sra

######Finally, change the file path of sequencing data of testing data in config_WES.yaml file.

