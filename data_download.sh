#!/bin/bash
# download and process databse file 
VEP_release="release-89"
echo "download reference file"
mkdir database && cd database
mkdir Fasta && cd Fasta
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz && gunzip Homo_sapiens_assembly38.fasta.gz
python ../../bin/reference_process.py
bwa index human.fasta
java -jar ../../software/picard.jar CreateSequenceDictionary R=human.fasta O=human.dict
samtools faidx human.fasta
python ../../bin/pyfasta_index.py 
rm Homo_sapiens_assembly38.fasta
cd ..

mkdir VCF_annotation && cd VCF_annotation
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_138.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_138.hg38.vcf.gz.tbi
cd ..

mkdir Protein && cd Protein
wget ftp://ftp.ensembl.org/pub/${VEP_release}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz && gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
mv Homo_sapiens.GRCh38.cdna.all.fa human.cdna.all.fa
wget ftp://ftp.ensembl.org/pub/${VEP_release}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz && gunzip Homo_sapiens.GRCh38.pep.all.fa.gz
mv Homo_sapiens.GRCh38.pep.all.fa human.pep.all.fa
makeblastdb -in human.pep.all.fa -dbtype prot -out peptide_database/peptide -parse_seqids
cd ..
cd ..






