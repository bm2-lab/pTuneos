import os,sys,time
import multiprocessing
import shutil 
import subprocess
import pandas as pd
import math
from pyper import *
import numpy as np
from sklearn import preprocessing
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from sklearn.semi_supervised import label_propagation
import itertools
from scipy import linalg
import matplotlib as mpl
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from imblearn.metrics import classification_report_imbalanced
from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from math import log, exp
from scipy import interp
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
import xgboost as xgb
from xgboost.sklearn import XGBClassifier
from sklearn import metrics
from sklearn.model_selection import GridSearchCV
import matplotlib.pylab as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier  
from imblearn.over_sampling import SMOTE
from collections import Counter
from sklearn.model_selection import cross_val_score
from sklearn.externals import joblib

a=26
k=4.86936
M=1. #default concentration of mutant peptides
W=1. #default concentration of wildtype peptides

WEPS=0.0003
HYDROPHOBIC_RESIDUES="AILMFWYV"
WEIRD_RESIDUES="CGP"
def get_iedb_seq(iedb_file):
	iedb_seq=[]
	for line in open(iedb_file):
		if line.startswith(">"):
			continue
		else:
			iedb_seq.append(line.strip())
	return iedb_seq
def read_trimmomatic(raw_fastq_path_first,raw_fastq_path_second,trimmomatic_path,adapter_path,fastq_prefix,logfile_fold,fastq_type,CPU):
	cmd_trimmomatic="java -jar " + trimmomatic_path + " PE -phred33 -threads " + str(CPU) + ' ' + raw_fastq_path_first + ' ' + raw_fastq_path_second  + ' -baseout ' + fastq_prefix + " ILLUMINACLIP:" + adapter_path + ':2:30:10' + ' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 > ' + logfile_fold + '/' + fastq_type + '_trimmomatic.log' + ' 2>&1'
	#print cmd_trimmomatic
	os.system(cmd_trimmomatic)


def hlatyping(raw_fastq_path_first,raw_fastq_path_second,opitype_fold,opitype_out_fold,opitype_ext,prefix,logfile_fold):
	cmd_hla = 'python ' + opitype_fold + ' -i ' + raw_fastq_path_first + ' ' + raw_fastq_path_second + ' --dna -o ' + opitype_out_fold + ' > ' + logfile_fold + '/' +  'hla_typing.log' + ' 2>&1'
	#print cmd_hla
	os.system(cmd_hla)
	result_dir=os.listdir(opitype_out_fold)
	#print result_dir[0]
	hla_result_path=opitype_out_fold+'/'+result_dir[0]+'/'+result_dir[0]+'_result.tsv'
	#print hla_result_path
	cmd_hla_ext = 'python ' + opitype_ext + ' -i ' + hla_result_path + ' -o ' + opitype_out_fold + ' -s ' + prefix
	#print cmd_hla_ext
	os.system(cmd_hla_ext)
	print 'hla type process done.'



def mapping_qc_gatk_preprocess(fastq_1_path,fastq_2_path,fastq_type,CPU,BWA_INDEX,alignment_out_fold,prefix,REFERENCE,bwa_path,samtools_path,java_picard_path,GATK_path,dbsnp138,OneKG,mills,logfile_fold,bamstat_out_fold):
	cmd_bwa=bwa_path + ' mem -t '+ str(CPU) + ' ' + BWA_INDEX + ' ' + fastq_1_path + ' ' +fastq_2_path + ' > ' + alignment_out_fold+'/'+'tmp_'+ prefix +'_'+fastq_type+'.sam'#+ logfile_fold + '/' + fastq_type + '_bwa.log' + ' 2>&1'
	cmd_samtools_1=samtools_path + ' view -bhS -@ '+ str(CPU) + ' ' + alignment_out_fold+'/'+'tmp_'+ prefix +'_'+fastq_type+'.sam' + ' -o ' + alignment_out_fold+'/'+'tmp_'+ prefix +'_'+fastq_type+'.bam > ' + logfile_fold + '/' + fastq_type + '_samtools_1.log' + ' 2>&1'
	cmd_samtools_sort=samtools_path + ' sort -@ ' + str(CPU) + ' -m 2G ' + alignment_out_fold+'/'+'tmp_'+ prefix +'_'+fastq_type+'.bam' + ' ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type + ' > ' + logfile_fold + '/' + fastq_type + '_samtools_sort.log' + ' 2>&1'
	cmd_samtools_index_1=samtools_path + ' index ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'.bam' + ' > ' + logfile_fold + '/' + fastq_type + '_samtools_index.log' + ' 2>&1'
	cmd_picard="java -Xmx4G -jar " + java_picard_path + ' MarkDuplicates INPUT=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'.bam' + ' OUTPUT=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_filter.bam' + ' METRICS_FILE=' + alignment_out_fold+'/'+prefix + '_'+fastq_type+'_dup_qc.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT > ' + logfile_fold + '/' + fastq_type + '_markdup.log' + ' 2>&1'
	cmd_samtools_index_2=samtools_path + ' index ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_filter.bam'
	cmd_add_readgroup="java -Xmx4G -jar " + java_picard_path + ' AddOrReplaceReadGroups I=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_filter.bam' + ' O=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_filter_add.bam' + ' SO=coordinate VALIDATION_STRINGENCY=SILENT RGID=' + fastq_type +  ' RGLB=' + fastq_type + ' RGPL=illumina RGSM='+ fastq_type + ' RGPU=NextSeq > ' + logfile_fold + '/' + fastq_type + '_addreadgroup.log' + ' 2>&1'
	cmd_buildbamindex="java -Xmx4G -jar " + java_picard_path + ' BuildBamIndex I=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_filter_add.bam' + ' O=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_filter_add.bam.bai' + ' VALIDATION_STRINGENCY=SILENT > ' + logfile_fold + '/' + fastq_type + '_buildindex.log' + ' 2>&1'
	cmd_BaseRecalibrator="java -Xmx4G -jar " + GATK_path + ' -T BaseRecalibrator -nct ' + str(CPU) + ' -R ' + REFERENCE + ' -I ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_filter_add.bam' + ' -knownSites ' + OneKG + ' -knownSites ' + mills + ' -knownSites ' + dbsnp138 + ' -o ' + alignment_out_fold + '/' + prefix + '_'+fastq_type + '.table > ' + logfile_fold + '/' + fastq_type + '_BaseRecalibrator.log' + ' 2>&1'
	cmd_PrintReads="java -Xmx4G -jar " + GATK_path + ' -T PrintReads -nct 8 -dt NONE -R ' + REFERENCE + ' -I ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_filter_add.bam' + ' -BQSR ' + alignment_out_fold + '/' + prefix + '_'+fastq_type + '.table' + ' -o ' + alignment_out_fold + '/' + prefix + '_'+fastq_type + '_recal.bam > ' + logfile_fold + '/' + fastq_type + '_PrintRead.log' + ' 2>&1'
	#print cmd_bwa
	os.system(cmd_bwa)
	#print cmd_samtools_1
	os.system(cmd_samtools_1)
	#print cmd_samtools_sort
	os.system(cmd_samtools_sort)
	#print cmd_samtools_index_1
	os.system(cmd_samtools_index_1)
	#print cmd_picard
	os.system(cmd_picard)
	#print cmd_samtools_index_2
	os.system(cmd_samtools_index_2)
	#print cmd_add_readgroup
	os.system(cmd_add_readgroup)
	#print cmd_buildbamindex
	os.system(cmd_buildbamindex)
	#print cmd_BaseRecalibrator
	os.system(cmd_BaseRecalibrator)
	#print cmd_PrintReads
	os.system(cmd_PrintReads)

def GATK_mutect2(GATK_path,REFERENCE,alignment_out_fold,prefix,CPU,dbsnp138,somatic_out_fold,vcftools_path,vep_path,vep_cache_path,netmhc_out_path,tumor_depth_cutoff,tumor_vaf_cutoff,normal_vaf_cutoff,itunes_bin_path,human_peptide_path,logfile_fold):
	cmd_GATK="java -Xmx4G -jar " + GATK_path + ' -T MuTect2 -nct ' + str(CPU) + ' -R ' + REFERENCE  + ' -I:tumor ' + alignment_out_fold + '/' + prefix + '_'+ 'tumor_recal.bam ' + '-I:normal ' + alignment_out_fold + '/' + prefix + '_'+ 'normal_recal.bam ' + '--dbsnp ' + dbsnp138 + ' -o ' + somatic_out_fold + '/' + prefix + '_'+ 'mutect2.vcf > ' + logfile_fold + '/' + prefix + '_mutect2.log' + ' 2>&1'
	print cmd_GATK
	os.system(cmd_GATK)
	cmd_mutation_filter='grep ' + "\'^#\|chr[1-9]\{0,1\}[0-9XY]\\{0,1\\}\\b\'" + ' ' + somatic_out_fold + '/' + prefix + '_'+ 'mutect2.vcf' + ' > ' + somatic_out_fold + '/' + prefix + '_' + 'mutect2_filter.vcf'
	#print cmd_mutation_filter
	os.system(cmd_mutation_filter)
	cmd_vcftools_snv=vcftools_path + " --vcf " + somatic_out_fold + '/' + prefix + '_'+ 'mutect2_filter.vcf' + " --remove-filtered-all --remove-indels --recode --recode-INFO-all --out " + somatic_out_fold + '/' + prefix + '_'+ 'SNVs_only'
	cmd_vcftools_indel=vcftools_path + " --vcf " + somatic_out_fold + '/' + prefix + '_'+ 'mutect2_filter.vcf' + " --remove-filtered-all --keep-only-indels --recode --recode-INFO-all --out " + somatic_out_fold + '/' + prefix + '_'+ 'INDELs_only'
	#print cmd_vcftools_snv
	#print cmd_vcftools_indel
	os.system(cmd_vcftools_snv)
	os.system(cmd_vcftools_indel)
	cmd_snv_filter="python " + itunes_bin_path + "/snv_filter.py -i " + somatic_out_fold + '/' + prefix + '_'+ 'SNVs_only.recode.vcf' + " -d " + str(tumor_depth_cutoff) + " -v " + str(tumor_vaf_cutoff) + " -n " + str(normal_vaf_cutoff) + " -o " + somatic_out_fold + " -s " + prefix
	#print cmd_snv_filter
	os.system(cmd_snv_filter)
	cmd_vep=vep_path + " -i " + somatic_out_fold + '/' + prefix + '_'+ 'SNVs_filter.vcf' + " --cache --dir " + vep_cache_path + " --dir_cache " + vep_cache_path + " --force_overwrite --canonical --symbol -o STDOUT --offline | filter_vep --ontology --filter \"CANONICAL is YES and Consequence is missense_variant\" -o " + somatic_out_fold + '/' + prefix + '_'+ 'snv_vep_ann.txt' + " --force_overwrite"
	#print cmd_vep
	os.system(cmd_vep)
	cmd_vep_snv_all=vep_path + " -i " + somatic_out_fold + '/' + prefix + '_'+ 'SNVs_filter.vcf' + " --cache --dir " + vep_cache_path + " --dir_cache " + vep_cache_path + " --force_overwrite --canonical --symbol -o STDOUT --offline | filter_vep --ontology --filter \"Consequence is missense_variant\" -o " + somatic_out_fold + '/' + prefix + '_'+ 'snv_vep_ann_all.txt' + " --force_overwrite"
	#print cmd_vep_snv_all
	os.system(cmd_vep_snv_all)
	cmd_vep_indel=vep_path + " -i " + somatic_out_fold + '/' + prefix + '_'+ 'INDELs_only.recode.vcf' + " --cache --dir " + vep_cache_path + " --dir_cache " + vep_cache_path + " --force_overwrite --canonical --symbol -o STDOUT --offline | filter_vep --ontology --filter \"Consequence is missense_variant\" -o " + somatic_out_fold + '/' + prefix + '_'+ 'mutect_indel_vep_ann.txt' + " --force_overwrite"
	#print cmd_vep_indel
	os.system(cmd_vep_indel)
	cmd_snv="python " + itunes_bin_path + "/snv2fasta.py -i " + somatic_out_fold + '/' + prefix + '_'+ 'snv_vep_ann.txt' + " -o " + netmhc_out_path + " -s " + prefix + " -p " + human_peptide_path
	#print cmd_snv
	os.system(cmd_snv)



def netMHCpan(fasta_file,hla_str,netmhc_out_file,out_dir,split_num,netMHCpan_path,tmp_dir,peptide_length):
	str_proc=r'''
#set -x
input_fasta=%s
hla_str=%s
netmhc_out=%s
out_dir=%s
split_num=%s
netMHCpan=%s
tmp=%s
pep_len=%s
if [ -d ${out_dir}/${tmp} ];then
	rm -rf ${out_dir}/${tmp}
	mkdir -p ${out_dir}/${tmp}
else
	mkdir -p ${out_dir}/${tmp}
fi
if [ -f ${netmhc_out} ];then
	rm ${netmhc_out}
fi
split -l ${split_num} ${input_fasta} ${out_dir}/${tmp}/
filelist=`ls ${out_dir}/${tmp}/`
arr1=(${filelist})
OLD_IFS="$IFS" 
IFS=","
arr2=(${hla_str})
IFS="$OLD_IFS" 
for s in ${arr2[@]}
do
{
	for file_l in ${arr1[@]}
	do
	{
		$netMHCpan -a $s -f ${out_dir}/${tmp}/${file_l} -l ${pep_len} -BA > ${out_dir}/${tmp}/${s}_${file_l}_tmp_netmhc.txt
	} &
	done
	wait
}
done
for file_l in ${arr1[@]}
do
{
	rm ${out_dir}/${tmp}/${file_l}
}
done
filelist1=`ls ${out_dir}/${tmp}/`
for file_r in $filelist1
do
{
	cat ${out_dir}/${tmp}/${file_r} >> ${netmhc_out}
	rm ${out_dir}/${tmp}/${file_r}	
}
done
rm -rf 	${out_dir}/${tmp}
#set +x
'''%(fasta_file,hla_str,netmhc_out_file,out_dir,split_num,netMHCpan_path,tmp_dir,peptide_length)
	subprocess.call(str_proc, shell=True, executable='/bin/bash')
def varscan_somatic_caling_drift(somatic_mutation_fold,alignment_out_fold,PREFIX,REFERENCE,vep_cache,samtools_path,varscan_path,vep_path,netmhc_out_fold,logfile_fold):
	str_proc = r'''
set -e
somat_f=%s
alignment_fold=%s
PREFIX=%s
REFERENCE=%s
vep_cache=%s
netmhc_out=%s
samtools=%s
varscan=%s
vep=%s
logfile_fold=%s
if [ ! -d ${somat_f} ];then
	mkdir ${somat_f}	
fi
if [ ! -d ${netmhc_out} ];then
	mkdir ${netmhc_out}	
fi
#rm -rf ${somat_f}/*
cd ${somat_f}
mkfifo ${PREFIX}_normal.fifo
mkfifo ${PREFIX}_tumor.fifo
$samtools mpileup -f ${REFERENCE} -q 5 -Q 20 -L 10000 -d 10000 ${alignment_fold}/${PREFIX}_normal_recal.bam  > ${PREFIX}_normal.fifo &
$samtools mpileup -f ${REFERENCE} -q 5 -Q 20 -L 10000 -d 10000 ${alignment_fold}/${PREFIX}_tumor_recal.bam  > ${PREFIX}_tumor.fifo &
java -jar $varscan somatic ${PREFIX}_normal.fifo ${PREFIX}_tumor.fifo ${PREFIX} > ${logfile_fold}/${PREFIX}_somatic.log 2>&1
java -jar $varscan processSomatic ${PREFIX}.snp
grep '^chrom\|chr[1-9]\{0,1\}[0-9XY]\{0,1\}\b' ${PREFIX}.snp > ${PREFIX}_filter.snp
rm ${PREFIX}_normal.fifo ${PREFIX}_tumor.fifo
cd ..
'''%(somatic_mutation_fold,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,samtools_path,varscan_path,vep_path,logfile_fold)
	#print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')


def varscan_neo(snv_fasta_file,hla_str,driver_gene_path,snv_netmhc_out_file,netmhc_out_fold,split_num,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_fold,netMHCpan_path,itunes_bin_path,peptide_length):
	netMHCpan(snv_fasta_file,hla_str,snv_netmhc_out_file,netmhc_out_fold,split_num,netMHCpan_path,'tmp_snv',peptide_length)
	str_proc1=r'''
PREFIX=%s
netmhc_out=%s
Exp_file=%s
Binding_Aff_Fc_Cutoff=%d
Binding_Aff_Cutoff=%d
Fpkm_Cutoff=%d
hla_str=%s
driver_gene_path=%s
netctl_fold=%s
itunes_bin_path=%s
python ${itunes_bin_path}/sm_netMHC_result_parse.py -i ${netmhc_out}/${PREFIX}_snv_netmhc.tsv -g ${netmhc_out}/${PREFIX}_snv.fasta -o ${netmhc_out} -s ${PREFIX}_snv -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
python ${itunes_bin_path}/netCTLPAN.py -i ${netmhc_out}/${PREFIX}_snv_final_neo_candidate.tsv -d ${driver_gene_path} -o ${netctl_fold} -s ${PREFIX}_snv
'''%(prefix,netmhc_out_fold,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,hla_str,driver_gene_path,netctl_fold,itunes_bin_path)
	#print str_proc1
	subprocess.call(str_proc1, shell=True, executable='/bin/bash')

def indel_calling_drift(strelka_out_fold,strelka_path,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,CPU,vep_path,itunes_bin_path):
	str_proc2=r'''
set -e
strelka_fold=%s
strelka_path=%s
alignment_fold=%s
PREFIX=%s
REFERENCE=%s
vep_cache=%s
netmhc_out=%s
cpu=%s
vep=%s
itunes_bin_path=%s
if [ -d ${strelka_fold} ];then
	rm -rf ${strelka_fold}
fi
python ${strelka_path}/configureStrelkaSomaticWorkflow.py --tumorBam=${alignment_fold}/${PREFIX}_tumor_recal.bam --normalBam=${alignment_fold}/${PREFIX}_normal_recal.bam --referenceFasta=${REFERENCE} --config=${strelka_path}/configureStrelkaSomaticWorkflow.py.ini --runDir=${strelka_fold} --exome
python ${strelka_fold}/runWorkflow.py -m local -j $cpu -q ${PREFIX}_strelka -g 8 --quiet
gunzip ${strelka_fold}/results/variants/somatic.indels.vcf.gz
$vep -i ${strelka_fold}/results/variants/somatic.indels.vcf --cache --dir ${vep_cache} --dir_cache ${vep_cache} --force_overwrite --canonical --symbol -o STDOUT --offline | filter_vep --ontology --filter "CANONICAL is YES and Consequence is coding_sequence_variant" -o ${strelka_fold}/${PREFIX}_strelka_indel_vep_ann.txt --force_overwrite
python ${itunes_bin_path}/varscandel2fasta.py -i ${strelka_fold}/${PREFIX}_strelka_indel_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}_strelka
python ${itunes_bin_path}/varscanins2fasta.py -i ${strelka_fold}/${PREFIX}_strelka_indel_vep_ann.txt  -o ${netmhc_out} -s ${PREFIX}_strelka
'''%(strelka_out_fold,strelka_path,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,CPU,vep_path,itunes_bin_path)
	#print str_proc2
	subprocess.call(str_proc2, shell=True, executable='/bin/bash')
	
def indel_neo(somatic_mutation_fold,PREFIX,vep_cache,netmhc_out_fold,vep_path,indel_fasta_file,hla_str,driver_gene_path,indel_netmhc_out_file,split_num,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_fold,netMHCpan_path,itunes_bin_path,peptide_length):
	str_proc1='''
somatic_mutation=%s
PREFIX=%s
vep_cache=%s
netmhc_out=%s
vep=%s
itunes_bin_path=%s
python ${itunes_bin_path}/varscan_indel_preprocess.py -i ${somatic_mutation}/${PREFIX}.indel -o ${somatic_mutation} -s ${PREFIX}
$vep -i ${somatic_mutation}/${PREFIX}_varscan_indel.vcf --cache --dir $vep_cache --dir_cache $vep_cache --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is coding_sequence_variant" -o ${somatic_mutation}/${PREFIX}_varscan_indel_vep_ann.txt --force_overwrite
python ${itunes_bin_path}/varscandel2fasta.py -i ${somatic_mutation}/${PREFIX}_varscan_indel_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}_varscan
python ${itunes_bin_path}/varscanins2fasta.py -i ${somatic_mutation}/${PREFIX}_varscan_indel_vep_ann.txt  -o ${netmhc_out} -s ${PREFIX}_varscan
'''%(somatic_mutation_fold,PREFIX,vep_cache,netmhc_out_fold,vep_path,itunes_bin_path)
	#print str_proc1
	subprocess.call(str_proc1, shell=True, executable='/bin/bash')	
	str_proc3=r'''
PREFIX=%s
netmhc_out=%s
cat ${netmhc_out}/${PREFIX}_strelka_del.fasta > ${netmhc_out}/${PREFIX}_indel.fasta
cat ${netmhc_out}/${PREFIX}_strelka_ins.fasta >> ${netmhc_out}/${PREFIX}_indel.fasta
cat ${netmhc_out}/${PREFIX}_varscan_del.fasta >> ${netmhc_out}/${PREFIX}_indel.fasta
cat ${netmhc_out}/${PREFIX}_varscan_ins.fasta >> ${netmhc_out}/${PREFIX}_indel.fasta
'''%(PREFIX,netmhc_out_fold)
	subprocess.call(str_proc3, shell=True, executable='/bin/bash')
	netMHCpan(indel_fasta_file,hla_str,indel_netmhc_out_file,netmhc_out_fold,split_num,netMHCpan_path,'tmp_indel',peptide_length)
	str_proc4=r'''
set -e
PREFIX=%s
netmhc_out=%s
Exp_file=%s
Binding_Aff_Fc_Cutoff=%s
Binding_Aff_Cutoff=%s
Fpkm_Cutoff=%s
hla_str=%s
driver_gene_path=%s
netctl_fold=%s
itunes_bin_path=%s
python ${itunes_bin_path}/sm_netMHC_result_parse.py -i ${netmhc_out}/${PREFIX}_indel_netmhc.tsv -g ${netmhc_out}/${PREFIX}_indel.fasta -o ${netmhc_out} -s ${PREFIX}_indel -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
python ${itunes_bin_path}/netCTLPAN.py -i ${netmhc_out}/${PREFIX}_indel_final_neo_candidate.tsv -d ${driver_gene_path} -o ${netctl_fold} -s ${PREFIX}_indel
'''%(PREFIX,netmhc_out_fold,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,hla_str,driver_gene_path,netctl_fold,itunes_bin_path)
	subprocess.call(str_proc4, shell=True, executable='/bin/bash')
	
	
def varscan_copynumber_calling(varscan_copynumber_fold,prefix,alignment_out_fold,REFERENCE,samtools_path,varscan_path,logfile_fold):
	str_proc=r'''
copynumber_profile=%s
PREFIX=%s
alignments=%s
REFERENCE=%s
samtools=%s
varscan=%s
logfile_fold=%s
if [ ! -d ${copynumber_profile} ];then
	mkdir ${copynumber_profile}	
fi
rm -rf ${copynumber_profile}/*
cd ${copynumber_profile}
mkfifo ${PREFIX}_normal.fifo
mkfifo ${PREFIX}_tumor.fifo
$samtools mpileup -f ${REFERENCE} -q 5 -Q 20 -L 10000 -d 10000 ${alignments}/${PREFIX}_normal_recal.bam > ${PREFIX}_normal.fifo &
$samtools mpileup -f ${REFERENCE} -q 5 -Q 20 -L 10000 -d 10000 ${alignments}/${PREFIX}_tumor_recal.bam > ${PREFIX}_tumor.fifo &
java -jar $varscan copynumber ${PREFIX}_normal.fifo ${PREFIX}_tumor.fifo ${PREFIX} > ${logfile_fold}/${PREFIX}_copynumber.log 2>&1
grep '^chrom\|chr[1-9]\{0,1\}[0-9XY]\{0,1\}\b' ${PREFIX}.copynumber > ${PREFIX}_filter.copynumber
rm ${PREFIX}_normal.fifo ${PREFIX}_tumor.fifo
cd ..
'''%(varscan_copynumber_fold,prefix,alignment_out_fold,REFERENCE,samtools_path,varscan_path,logfile_fold)
	subprocess.call(str_proc, shell=True, executable='/bin/bash')

def pyclone_annotation(somatic_mutation_fold,varscan_copynumber_fold,prefix,pyclone_fold,netctl_fold,pyclone_path,logfile_fold,itunes_bin_path):
	str_proc=r'''
somatic_mutation=%s
copynumber_profile=%s
PREFIX=%s
pyclone=%s
netctl=%s
Pyclone=%s
logfile_fold=%s
itunes_bin_path=%s
Rscript ${itunes_bin_path}/sequenza_test.R ${somatic_mutation}/${PREFIX}_filter.snp ${copynumber_profile}/${PREFIX}_filter.copynumber ${copynumber_profile}/ ${PREFIX} > ${logfile_fold}/${PREFIX}_sequenza.log 2>&1
python ${itunes_bin_path}/pyclone_input.py -n ${netctl}/${PREFIX}_snv_netctl_concact.tsv -i ${somatic_mutation}/${PREFIX}_snv_vep_ann_all.txt -s ${somatic_mutation}/${PREFIX}_SNVs_only.recode.vcf -c ${copynumber_profile}/${PREFIX}_seg_copynumber.txt -o ${pyclone} -S ${PREFIX}
TUMOR_CONTENT=`cat ${copynumber_profile}/${PREFIX}_cellularity.txt`
$Pyclone setup_analysis --in_files ${pyclone}/${PREFIX}_pyclone_input.tsv --tumour_contents $TUMOR_CONTENT --prior major_copy_number --working_dir ${pyclone}
$Pyclone run_analysis --config_file ${pyclone}/config.yaml > ${logfile_fold}/${PREFIX}_pyclone.log 2>&1
$Pyclone build_table --config_file ${pyclone}/config.yaml --out_file ${pyclone}/loci.tsv --table_type loci
python ${itunes_bin_path}/neo_pyclone_annotation.py -n ${netctl}/${PREFIX}_snv_netctl_concact.tsv -i ${somatic_mutation}/${PREFIX}_snv_vep_ann_all.txt -s ${pyclone}/loci.tsv -o ${netctl} -S ${PREFIX}
'''%(somatic_mutation_fold,varscan_copynumber_fold,prefix,pyclone_fold,netctl_fold,pyclone_path,logfile_fold,itunes_bin_path)
	#print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')


def kallisto_expression(raw_fastq_path_first,raw_fastq_path_second,kallisto_path,kallisto_out_fold,prefix,kallisto_cdna_path,logfile_fold):
	cdna_path_dir = os.path.dirname(kallisto_cdna_path)
	cnd_file_prefix = os.path.splitext(os.path.basename(kallisto_cdna_path))[0]
	kallisto_index_path = cdna_path_dir + '/' + cnd_file_prefix + '.idx'
	cmd_kallisto_index = kallisto_path + " index -i " + kallisto_index_path + ' ' + kallisto_cdna_path + ' > ' +  logfile_fold + '/' + prefix + '_kallisto_index.log' + ' 2>&1'
	cmd_kallisto_quant = kallisto_path + " quant -i " + kallisto_index_path + " -o " + kallisto_out_fold + " " + raw_fastq_path_first + " " + raw_fastq_path_second + ' > ' +  logfile_fold + '/' + prefix + '_kallisto.log' + ' 2>&1'
	#print cmd_kallisto_index
	os.system(cmd_kallisto_index)
	#print cmd_kallisto_quant
	os.system(cmd_kallisto_quant)


def hydro_vector(pep):
	hydro_score={"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,"K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,"T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3}
	hydrophobicity_vector=[]
	pep_list=list(pep)
	pep_len=len(pep_list)
	for pep in pep_list:
		hydrophobicity_vector.append(hydro_score[pep.upper()])
	return hydrophobicity_vector



def logSum(v):
	ma=max(v)
	return log(sum(map(lambda x: exp(x-ma),v)))+ma

def aligner(seq1,seq2):
	matrix = matlist.blosum62
	gap_open = -11
	gap_extend = -1
	aln = pairwise2.align.localds(seq1.upper(), seq2.upper(), matrix,gap_open,gap_extend)
	#print aln
	return aln

def calculate_R(neo_seq,iedb_seq):
	align_score=[]
	#i=0
	for seq in iedb_seq:
		aln_score=aligner(neo_seq,seq)
		#i=i+1
		#print i
		#print aln_score
		if aln_score!=[]:
			localds_core=max([line[2] for line in aln_score])
			align_score.append(localds_core)
	#print align_score
	#print k,a
	bindingEnergies=map(lambda x: -k*(a-x),align_score)
	#print bindingEnergies
	lZk=logSum(bindingEnergies+[0])
	lGb=logSum(bindingEnergies)
	R=exp(lGb-lZk)
	return R

############similarity############

def cal_similarity_per(mut_seq,normal_seq):
	score_pair=aligner(mut_seq,normal_seq)[0][2]
	score_self=aligner(mut_seq,mut_seq)[0][2]
	per_similarity=score_pair/score_self
	return per_similarity


def get_homolog_info(mut_seq,hla_type,blastp_tmp_file,blastp_out_tmp_file,netMHCpan_pep_tmp_file,netMHCpan_ml_out_tmp_file,blast_db_path):
	hla_type_in=hla_type.replace('*','')
	blastp_fasta_line='>'+'\n'+mut_seq+'\n'
	#print blastp_fasta_line
	pep_len=len(mut_seq)
	f=open(blastp_tmp_file,'w')
	f.write(blastp_fasta_line)
	f.close()
	str_blastp_pro='blastp -query ' + blastp_tmp_file + ' -db ' + blast_db_path + ' -out ' + blastp_out_tmp_file + ' -evalue 200000 -comp_based_stats 0'
	#print str_blastp_pro
	subprocess.call(str_blastp_pro,shell = True,executable = '/bin/bash')
	for line in open(blastp_out_tmp_file):
		if line.startswith('Sbjct'): 
			human_pep_record=line.strip().split(' ')
			human_pep = [i for i in human_pep_record if i!=''][2]
			#print len(human_pep)
			if len(human_pep)==pep_len:
				human_homolog_pep=(human_pep)
				break
			else:
				continue
		else:
			continue
	if len(human_pep)<pep_len:
		human_homolog_pep=mut_seq
	#print human_homolog_pep
	f=open(netMHCpan_pep_tmp_file,'w')
	f.write(human_homolog_pep+'\n')
	f.close()
	str_netMHCpan_ml_pro='netMHCpan -p ' + netMHCpan_pep_tmp_file + ' -a ' + hla_type_in + ' > ' + netMHCpan_ml_out_tmp_file
	#print str_netMHCpan_ml_pro
	subprocess.call(str_netMHCpan_ml_pro,shell = True,executable = '/bin/bash')
	for line in open(netMHCpan_ml_out_tmp_file):
		if not line.startswith('    '):
			continue
		else:
			record=line.strip().split(' ')
			ml_record = [i for i in record if i!='']
			human_homolog_pep_el=ml_record[12]
	return human_homolog_pep,human_homolog_pep_el

def get_EL_info(seq,hla_type,netMHCpan_pep_tmp_file,netMHCpan_ml_out_tmp_file):
	hla_type_in=hla_type.replace('*','')
	f=open(netMHCpan_pep_tmp_file,'w')
	f.write(seq+'\n')
	f.close()
	str_netMHCpan_ml_pro='netMHCpan -p ' + netMHCpan_pep_tmp_file + ' -a ' + hla_type_in + ' > ' + netMHCpan_ml_out_tmp_file
	#print str_netMHCpan_ml_pro
	subprocess.call(str_netMHCpan_ml_pro,shell = True,executable = '/bin/bash')
	pep_el_rank=[]
	for line in open(netMHCpan_ml_out_tmp_file):
		if not line.startswith('    '):
			continue
		else:
			record=line.strip().split(' ')
			ml_record = [i for i in record if i!='']
			pep_el_rank=ml_record[12]
	return pep_el_rank
#read in data 
def InVivoModelAndScoreSNV(neo_file,cf_hy_model_9,cf_hy_model_10,cf_hy_model_11,RF_model,neo_model_file,blastp_tmp_file,blastp_out_tmp_file,netMHCpan_pep_tmp_file,netMHCpan_ml_out_tmp_file,iedb_file,blast_db_path):
	iedb_seq=get_iedb_seq(iedb_file)
	hy_xgb_9=joblib.load(cf_hy_model_9)
	hy_xgb_10=joblib.load(cf_hy_model_10)
	hy_xgb_11=joblib.load(cf_hy_model_11)
	#print hy_xgb_9
	data_neo=pd.read_table(neo_file,header=0,sep='\t')
	MT_peptide=data_neo.MT_pep
	HLA=data_neo.HLA_type
	WT_peptide=data_neo.WT_pep
	hydrophobicity_score=[]
	Recognition_score=[]
	Homolog_pep=[]
	Homolog_EL=[]
	MT_peptide_EL=[]
	WT_peptide_EL=[]
	for i in range(len(MT_peptide)):
		line=MT_peptide[i]
		H_p,H_E=get_homolog_info(line,HLA[i],blastp_tmp_file,blastp_out_tmp_file,netMHCpan_pep_tmp_file,netMHCpan_ml_out_tmp_file,blast_db_path)
		mt_el=get_EL_info(MT_peptide[i],HLA[i],netMHCpan_pep_tmp_file,netMHCpan_ml_out_tmp_file)
		wt_el=get_EL_info(WT_peptide[i],HLA[i],netMHCpan_pep_tmp_file,netMHCpan_ml_out_tmp_file)
		MT_peptide_EL.append(mt_el)
		WT_peptide_EL.append(wt_el)
		Homolog_pep.append(H_p)
		Homolog_EL.append(H_E)
		if len(line)==9:
			h_score=hy_xgb_9.predict_proba(np.array(hydro_vector(line)).reshape((1,9)))[:,1][0]
			hydrophobicity_score.append(h_score)
			R=calculate_R(line,iedb_seq)
			Recognition_score.append(R)
		elif len(line)==10:
			h_score=hy_xgb_10.predict_proba(np.array(hydro_vector(line)).reshape((1,10)))[:,1][0]
			hydrophobicity_score.append(h_score)
			R=calculate_R(line,iedb_seq)
			Recognition_score.append(R)
		elif len(line)==11:
			h_score=hy_xgb_11.predict_proba(np.array(hydro_vector(line)).reshape((1,11)))[:,1][0]
			hydrophobicity_score.append(h_score)
			R=calculate_R(line,iedb_seq)
			Recognition_score.append(R)
		else:
			print "The length of peptide is out of our considertion!!"
			hydrophobicity_score.append(0.5)
			R=calculate_R(line,iedb_seq)	
			Recognition_score.append(R)
	paired_similarity_score=[]
	homolog_similaity_score=[]
	#####paired similarity and homolog similarity########
	for M_P,N_P,H_P in zip(data_neo.MT_pep,data_neo.WT_pep,Homolog_pep):
		paired_s=cal_similarity_per(M_P,N_P)
		homolog_s=cal_similarity_per(M_P,H_P)
		paired_similarity_score.append(paired_s)
		homolog_similaity_score.append(homolog_s)
	self_sequence_similarity=[]
	for i in range(len(paired_similarity_score)):
		if paired_similarity_score[i] >= homolog_similaity_score[i]:
			sss=paired_similarity_score[i]
		else:
			sss=homolog_similaity_score[i]
		self_sequence_similarity.append(sss)
	data_neo["Homolog_pep"]=Homolog_pep
	data_neo["Homolog_Binding_EL"]=Homolog_EL
	data_neo["Recognition_score"]=Recognition_score
	data_neo["Hydrophobicity_score"]=hydrophobicity_score
	data_neo["Self_sequence_similarity"]=self_sequence_similarity
	data_neo["MT_Binding_EL"]=MT_peptide_EL
	data_neo["WT_Binding_EL"]=WT_peptide_EL
	df_neo=data_neo.loc[:,['Hydrophobicity_score','Recognition_score','Self_sequence_similarity','MT_Binding_EL','WT_Binding_EL']]
	cf_RF=joblib.load(RF_model)
	dneo_predprob = cf_RF.predict_proba(df_neo.values)[:,1]
	data_neo["model_pro"]=dneo_predprob
	f_EL_rank_wt=lambda x:1-(1/(1+math.pow(math.e,5*(float(x)-2))))/2
	f_EL_rank_mt=lambda x:1/(1+math.pow(math.e,5*(float(x)-2)))
	EL_mt_rank_score=data_neo.MT_Binding_EL.apply(f_EL_rank_mt)
	EL_wt_rank_score=data_neo.WT_Binding_EL.apply(f_EL_rank_wt)
	k=1
	f_TPM=lambda x:math.tanh(x/k)
	allele_frequency_score=data_neo.variant_allele_frequency
	netchop_score=data_neo.combined_prediction_score
	cellular_prevalence_score=data_neo.cellular_prevalence
	tpm_score=data_neo.tpm.apply(f_TPM)
	immuno_effect_score=[tpm_score[i]*allele_frequency_score[i]*netchop_score[i]*cellular_prevalence_score[i]*data_neo.Hydrophobicity_score[i]*data_neo.Recognition_score[i]*data_neo.Self_sequence_similarity[i]*EL_mt_rank_score[i]*EL_wt_rank_score[i] for i in range(len(data_neo.MT_Binding_EL))]
	data_neo["immuno_effect_score"]=immuno_effect_score
	data_neo_out_sort=data_neo.sort_values(['model_pro',"immuno_effect_score"],ascending=[0,0])
	del data_neo_out_sort["contain_X"]
	del data_neo_out_sort["target_id"]
	data_neo_out_sort.to_csv(neo_model_file,sep='\t',header=1,index=0)
		

def InVivoModelAndScoreINDEL(neo_file,cf_hy_model_9,cf_hy_model_10,cf_hy_model_11,RF_model,neo_model_file,blastp_tmp_file,blastp_out_tmp_file,netMHCpan_pep_tmp_file,netMHCpan_ml_out_tmp_file,iedb_file,blast_db_path):
	iedb_seq=get_iedb_seq(iedb_file)
	hy_xgb_9=joblib.load(cf_hy_model_9)
	hy_xgb_10=joblib.load(cf_hy_model_10)
	hy_xgb_11=joblib.load(cf_hy_model_11)
	data_neo=pd.read_table(neo_file,header=0,sep='\t')
	MT_peptide=data_neo.MT_pep
	HLA=data_neo.HLA_type
	WT_peptide=data_neo.WT_pep
	hydrophobicity_score=[]
	Recognition_score=[]
	Homolog_pep=[]
	Homolog_EL=[]
	MT_peptide_EL=[]
	WT_peptide_EL=[]
	for i in range(len(MT_peptide)):
		line=MT_peptide[i]
		H_p,H_E=get_homolog_info(line,HLA[i],blastp_tmp_file,blastp_out_tmp_file,netMHCpan_pep_tmp_file,netMHCpan_ml_out_tmp_file,blast_db_path)
		mt_el=get_EL_info(MT_peptide[i],HLA[i],netMHCpan_pep_tmp_file,netMHCpan_ml_out_tmp_file)
		wt_el=get_EL_info(WT_peptide[i],HLA[i],netMHCpan_pep_tmp_file,netMHCpan_ml_out_tmp_file)
		MT_peptide_EL.append(mt_el)
		WT_peptide_EL.append(wt_el)
		Homolog_pep.append(H_p)
		Homolog_EL.append(H_E)
		if len(line)==9:
			h_score=hy_xgb_9.predict_proba(np.array(hydro_vector(line)).reshape((1,9)))[:,1][0]
			hydrophobicity_score.append(h_score)
			R=calculate_R(line,iedb_seq)
			Recognition_score.append(R)
		elif len(line)==10:
			h_score=hy_xgb_10.predict_proba(np.array(hydro_vector(line)).reshape((1,10)))[:,1][0]
			hydrophobicity_score.append(h_score)
			R=calculate_R(line,iedb_seq)
			Recognition_score.append(R)
		elif len(line)==11:
			h_score=hy_xgb_11.predict_proba(np.array(hydro_vector(line)).reshape((1,11)))[:,1][0]
			hydrophobicity_score.append(h_score)
			R=calculate_R(line,iedb_seq)
			Recognition_score.append(R)
		else:
			print "The length of peptide is out of our considertion!!"
			hydrophobicity_score.append(0.5)
			R=calculate_R(line,iedb_seq)	
			Recognition_score.append(R)
	paired_similarity_score=[]
	homolog_similaity_score=[]
	#####paired similarity and homolog similarity########
	for M_P,N_P,H_P in zip(data_neo.MT_pep,data_neo.WT_pep,Homolog_pep):
		paired_s=cal_similarity_per(M_P,N_P)
		homolog_s=cal_similarity_per(M_P,H_P)
		paired_similarity_score.append(paired_s)
		homolog_similaity_score.append(homolog_s)
	self_sequence_similarity=[]
	for i in range(len(paired_similarity_score)):
		if paired_similarity_score[i] >= homolog_similaity_score[i]:
			sss=paired_similarity_score[i]
		else:
			sss=homolog_similaity_score[i]
		self_sequence_similarity.append(sss)
	data_neo["Homolog_pep"]=Homolog_pep
	data_neo["Homolog_Binding_EL"]=Homolog_EL
	data_neo["Recognition_score"]=Recognition_score
	data_neo["Hydrophobicity_score"]=hydrophobicity_score
	data_neo["Self_sequence_similarity"]=self_sequence_similarity
	data_neo["MT_Binding_EL"]=MT_peptide_EL
	data_neo["WT_Binding_EL"]=WT_peptide_EL
	df_neo=data_neo.loc[:,['Hydrophobicity_score','Recognition_score','Self_sequence_similarity','MT_Binding_EL','WT_Binding_EL']]
	cf_RF=joblib.load(RF_model)
	dneo_predprob = cf_RF.predict_proba(df_neo.values)[:,1]
	data_neo["model_pro"]=dneo_predprob
	f_EL_rank_wt=lambda x:1-(1/(1+math.pow(math.e,5*(float(x)-2))))/2
	f_EL_rank_mt=lambda x:1/(1+math.pow(math.e,5*(float(x)-2)))
	EL_mt_rank_score=data_neo.MT_Binding_EL.apply(f_EL_rank_mt)
	EL_wt_rank_score=data_neo.WT_Binding_EL.apply(f_EL_rank_wt)
	k=1
	f_TPM=lambda x:math.tanh(x/k)
	netchop_score=data_neo.combined_prediction_score
	tpm_score=data_neo.tpm.apply(f_TPM)
	immuno_effect_score=[tpm_score[i]*netchop_score[i]*data_neo.Hydrophobicity_score[i]*data_neo.Recognition_score[i]*data_neo.Self_sequence_similarity[i]*EL_mt_rank_score[i]*EL_wt_rank_score[i] for i in range(len(data_neo.MT_Binding_EL))]
	data_neo["immuno_effect_score"]=immuno_effect_score
	data_neo_out_sort=data_neo.sort_values(['model_pro',"immuno_effect_score"],ascending=[0,0])
	del data_neo_out_sort["contain_X"]
	del data_neo_out_sort["target_id"]
	data_neo_out_sort.to_csv(neo_model_file,sep='\t',header=1,index=0)
		