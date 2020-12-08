#!/bin/bash
# download and process databse file 
VEP_release="release-97"
echo "download reference file"
mkdir database && cd database
mkdir Fasta && cd Fasta
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.gz && gunzip ucsc.hg19.fasta.gz
mv ucsc.hg19.fasta human.fasta
bwa index human.fasta
java -jar ../../software/picard.jar CreateSequenceDictionary R=human.fasta O=human.dict
samtools faidx human.fasta
python ../../bin/pyfasta_index.py 
sequenza-utils gc_wiggle -w 50 --fasta human.fasta -o hg19.gc50Base.wig.gz
rm ucsc.hg19.fasta
cd ..

mkdir VCF_annotation && cd VCF_annotation
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
tabix Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
tabix 1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz
tabix dbsnp_138.hg19.vcf.gz
cd ..

mkdir Protein && cd Protein
wget ftp://ftp.ensembl.org/pub/grch37/${VEP_release}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz && Homo_sapiens.GRCh37.cdna.all.fa.gz
mv Homo_sapiens.GRCh37.cdna.all.fa human.cdna.all.fa
wget ftp://ftp.ensembl.org/pub/grch37/${VEP_release}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.pep.all.fa.gz && gunzip Homo_sapiens.GRCh37.pep.all.fa.gz
mv Homo_sapiens.GRCh37.pep.all.fa human.pep.all.fa
makeblastdb -in human.pep.all.fa -dbtype prot -out peptide_database/peptide -parse_seqids
cd ..
cd ..

mkdir vep_data && cd vep_data
wget ftp://ftp.ensembl.org/pub/grch37/${VEP_release}/variation/VEP/homo_sapiens_vep_97_GRCh37.tar.gz && tar xvzf homo_sapiens_vep_97_GRCh37.tar.gz












