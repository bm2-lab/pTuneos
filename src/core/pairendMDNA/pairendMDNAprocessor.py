import os,sys,time
import multiprocessing
import subprocess
import pandas as pd
import math
#from pyper import *
import numpy as np
from sklearn import preprocessing
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from sklearn.semi_supervised import label_propagation
import itertools
import matplotlib as mpl
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
#from imblearn.metrics import classification_report_imbalanced
from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from math import log, exp
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
import xgboost as xgb
from xgboost.sklearn import XGBClassifier
from sklearn import metrics
from sklearn.model_selection import GridSearchCV
import matplotlib.pylab as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier  
#from imblearn.over_sampling import SMOTE
from collections import Counter
from sklearn.model_selection import cross_val_score
from sklearn.externals import joblib
import warnings
warnings.filterwarnings("ignore")
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
	cmd_samtools_sort=samtools_path + ' sort -@ ' + '4 -m 2G ' + alignment_out_fold+'/'+'tmp_'+ prefix +'_'+fastq_type+'.bam' + ' ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type + ' > ' + logfile_fold + '/' + fastq_type + '_samtools_sort.log' + ' 2>&1'
	cmd_samtools_index_1=samtools_path + ' index ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'.bam' + ' > ' + logfile_fold + '/' + fastq_type + '_samtools_index.log' + ' 2>&1'
	cmd_picard="java -Xmx4G -jar " + java_picard_path + ' MarkDuplicates INPUT=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'.bam' + ' OUTPUT=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_filter.bam' + ' METRICS_FILE=' + alignment_out_fold+'/'+prefix + '_'+fastq_type+'_dup_qc.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT > ' + logfile_fold + '/' + fastq_type + '_markdup.log' + ' 2>&1'
	cmd_samtools_index_2=samtools_path + ' index ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_filter.bam'
	cmd_add_readgroup="java -Xmx4G -jar " + java_picard_path + ' AddOrReplaceReadGroups I=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_filter.bam' + ' O=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_filter_add.bam' + ' SO=coordinate VALIDATION_STRINGENCY=SILENT RGID=' + fastq_type +  ' RGLB=' + fastq_type + ' RGPL=illumina RGSM='+ fastq_type + ' RGPU=NextSeq > ' + logfile_fold + '/' + fastq_type + '_addreadgroup.log' + ' 2>&1'
	cmd_buildbamindex="java -Xmx4G -jar " + java_picard_path + ' BuildBamIndex I=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_filter_add.bam' + ' O=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_filter_add.bam.bai' + ' VALIDATION_STRINGENCY=SILENT > ' + logfile_fold + '/' + fastq_type + '_buildindex.log' + ' 2>&1'
	cmd_BaseRecalibrator="java -Xmx4G -jar " + GATK_path + ' -T BaseRecalibrator -nct ' + str(CPU) + ' -R ' + REFERENCE + ' -I ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_filter_add.bam' + ' -knownSites ' + OneKG + ' -knownSites ' + mills + ' -knownSites ' + dbsnp138 + ' -o ' + alignment_out_fold + '/' + prefix + '_'+fastq_type + '.table > ' + logfile_fold + '/' + fastq_type + '_BaseRecalibrator.log' + ' 2>&1'
	cmd_PrintReads="java -Xmx4G -jar " + GATK_path + ' -T PrintReads -nct 16 -dt NONE -R ' + REFERENCE + ' -I ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_filter_add.bam' + ' -BQSR ' + alignment_out_fold + '/' + prefix + '_'+fastq_type + '.table' + ' -o ' + alignment_out_fold + '/' + prefix + '_'+fastq_type + '_recal.bam > ' + logfile_fold + '/' + fastq_type + '_PrintRead.log' + ' 2>&1'
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

def GATK_mutect2(GATK_path,REFERENCE,alignment_out_fold,prefix,CPU,dbsnp138,somatic_out_fold,vcftools_path,vep_path,vep_cache_path,netmhc_out_path,tumor_depth_cutoff,tumor_vaf_cutoff,normal_vaf_cutoff,pTuneos_bin_path,human_peptide_path,logfile_fold):
	str_proc_gatk=r'''
set -e
GATK_path=%s
REFERENCE=%s
alignment=%s
logfile_fold=%s
somatic=%s
prefix=%s
dbsnp=%s
chr=(`seq 1 22` X Y)
for i in ${chr[@]}
do
{
	java -Xmx4G -jar ${GATK_path} -T MuTect2 -L chr${i} -R $REFERENCE -I:tumor ${alignment}/${prefix}_tumor_recal.bam -I:normal ${alignment}/${prefix}_normal_recal.bam --dbsnp $dbsnp -o ${alignment}/${i}.vcf > ${logfile_fold}/${prefix}_mutect2_${i}.log 2>&1
}&
done
wait
bcftools concat -o ${somatic}/${prefix}_mutect2.vcf \
${alignment}/1.vcf ${alignment}/2.vcf ${alignment}/3.vcf ${alignment}/4.vcf ${alignment}/5.vcf ${alignment}/6.vcf ${alignment}/7.vcf ${alignment}/8.vcf ${alignment}/9.vcf ${alignment}/10.vcf \
${alignment}/11.vcf ${alignment}/12.vcf ${alignment}/13.vcf ${alignment}/14.vcf ${alignment}/15.vcf ${alignment}/16.vcf ${alignment}/17.vcf ${alignment}/18.vcf \
${alignment}/19.vcf ${alignment}/20.vcf ${alignment}/21.vcf ${alignment}/22.vcf ${alignment}/X.vcf ${alignment}/Y.vcf
rm ${alignment}/1.vcf ${alignment}/2.vcf ${alignment}/3.vcf ${alignment}/4.vcf ${alignment}/5.vcf ${alignment}/6.vcf ${alignment}/7.vcf ${alignment}/8.vcf ${alignment}/9.vcf ${alignment}/10.vcf \
${alignment}/11.vcf ${alignment}/12.vcf ${alignment}/13.vcf ${alignment}/14.vcf ${alignment}/15.vcf ${alignment}/16.vcf ${alignment}/17.vcf ${alignment}/18.vcf \
${alignment}/19.vcf ${alignment}/20.vcf ${alignment}/21.vcf ${alignment}/22.vcf ${alignment}/X.vcf ${alignment}/Y.vcf
rm -rf {alignment}/*.idx
rm -rf {alignment}/*.vcf
	'''%(GATK_path,REFERENCE,alignment_out_fold,logfile_fold,somatic_out_fold,prefix,dbsnp138)
	print str_proc_gatk
	subprocess.call(str_proc_gatk, shell=True, executable='/bin/bash')
	cmd_mutation_filter='grep ' + "\'^#\|chr[1-9]\{0,1\}[0-9XY]\\{0,1\\}\\b\'" + ' ' + somatic_out_fold + '/' + prefix + '_'+ 'mutect2.vcf' + ' > ' + somatic_out_fold + '/' + prefix + '_' + 'mutect2_filter.vcf'
	#print cmd_mutation_filter
	os.system(cmd_mutation_filter)
	cmd_vcftools_pass=vcftools_path + " --vcf " + somatic_out_fold + '/' + prefix + '_' + 'mutect2_filter.vcf' + " --remove-filtered-all --recode --recode-INFO-all --out " + somatic_out_fold + '/' + prefix + '_'+ 'mutect2_pass > ' + logfile_fold + '/' + prefix + '_vcftools_snv.log 2>&1'
	os.system(cmd_vcftools_pass)
	cmd_filter="python " + pTuneos_bin_path + "/snv_filter.py -i " + somatic_out_fold + '/' + prefix + '_'+ 'mutect2_pass.recode.vcf' + " -d " + str(tumor_depth_cutoff) + " -v " + str(tumor_vaf_cutoff) + " -n " + str(normal_vaf_cutoff) + " -o " + somatic_out_fold + " -s " + prefix
	#print cmd_filter
	os.system(cmd_filter)
	cmd_vep_snv_all=vep_path + " -i " + somatic_out_fold + '/' + prefix + '_'+ 'filter.vcf' + " --cache --dir " + vep_cache_path + " --dir_cache " + vep_cache_path + " --force_overwrite --canonical --symbol -o STDOUT --offline | filter_vep --ontology --filter \"Consequence is coding_sequence_variant\" -o " + somatic_out_fold + '/' + prefix + '_'+ 'vep_ann_all.txt' + " --force_overwrite"
	#print cmd_vep_snv_all
	os.system(cmd_vep_snv_all)
	cmd_snv_peptide="python " + pTuneos_bin_path + "/snv2fasta.py -i " + somatic_out_fold + '/' + prefix + '_'+ 'vep_ann_all.txt' + ' -o ' + somatic_out_fold + ' -s ' + prefix + ' -p ' + human_peptide_path
	#print cmd_snv_peptide
	os.system(cmd_snv_peptide)
	cmd_insetion_peptide="python " + pTuneos_bin_path + "/varscanins2fasta.py -i " + somatic_out_fold + '/' + prefix + '_'+ 'vep_ann_all.txt' + ' -o ' + somatic_out_fold + ' -s ' + prefix + ' -r ' + REFERENCE +' -p ' + human_peptide_path
	#print cmd_insetion_peptide
	os.system(cmd_insetion_peptide)
	cmd_deletion_peptide="python " + pTuneos_bin_path + "/varscandel2fasta.py -i " + somatic_out_fold + '/' + prefix + '_'+ 'vep_ann_all.txt' + ' -o ' + somatic_out_fold + ' -s ' + prefix + ' -r ' + REFERENCE +' -p ' + human_peptide_path
	#print cmd_deletion_peptide
	os.system(cmd_deletion_peptide)
	cmd_cat_snv="cat " + somatic_out_fold + '/' + prefix + '_'+ 'snv.fasta > ' + netmhc_out_path + '/' +  prefix + '_'+ 'all.fasta'
	#print cmd_cat_snv
	os.system(cmd_cat_snv)
	cmd_cat_ins="cat " + somatic_out_fold + '/' + prefix + '_'+ 'ins.fasta >> ' + netmhc_out_path + '/' +  prefix + '_'+ 'all.fasta'
	#print cmd_cat_ins
	os.system(cmd_cat_ins)
	cmd_cat_del="cat " + somatic_out_fold + '/' + prefix + '_'+ 'del.fasta >> ' + netmhc_out_path + '/' +  prefix + '_'+ 'all.fasta'
	#print cmd_cat_del
	os.system(cmd_cat_del)



def netMHCpan(fasta_file,hla_str,netmhc_out_file,out_dir,split_num,netMHCpan_path,tmp_dir,peptide_length):
	str_proc=r'''
set -e
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
'''%(fasta_file,hla_str,netmhc_out_file,out_dir,split_num,netMHCpan_path,tmp_dir,peptide_length)
	subprocess.call(str_proc, shell=True, executable='/bin/bash')



def sequenza_cal(alignment_out_fold,sequenza_path,REFERENCE,gc_file_path,copynumber_fold,prefix,pTuneos_bin_path):
	str_proc=r'''
alignment_out_fold=%s
sequenza_path=%s
REFERENCE=%s
gc_file_path=%s
copynumber_fold=%s
PREFIX=%s
pTuneos_bin_path=%s
${sequenza_path} bam2seqz -n ${alignment_out_fold}/${PREFIX}_normal_recal.bam -t ${alignment_out_fold}/${PREFIX}_tumor_recal.bam --fasta ${REFERENCE} -gc ${gc_file_path} -o ${copynumber_fold}/${PREFIX}.out.seqz.gz
${sequenza_path} seqz_binning --seqz ${copynumber_fold}/${PREFIX}.out.seqz.gz -w 50 -o ${copynumber_fold}/${PREFIX}.small.seqz.gz
Rscript ${pTuneos_bin_path}/sequenza_process.R ${copynumber_fold}/${PREFIX}.small.seqz.gz ${copynumber_fold} ${PREFIX}
'''%(alignment_out_fold,sequenza_path,REFERENCE,gc_file_path,copynumber_fold,prefix,pTuneos_bin_path)
	#print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')






def neo_cal(all_fasta_file,hla_str,driver_gene_path,netmhc_out_file,netmhc_out_fold,split_num,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_out_fold,netMHCpan_path,peptide_length,pTuneos_bin_path,netchop_path):
	netMHCpan(all_fasta_file,hla_str,netmhc_out_file,netmhc_out_fold,split_num,netMHCpan_path,'tmp_neo',peptide_length)
	str_proc1=r'''
PREFIX=%s
netmhc_out=%s
Exp_file=%s
Binding_Aff_Fc_Cutoff=%d
Binding_Aff_Cutoff=%d
Fpkm_Cutoff=%d
hla_str=%s
driver_gene_path=%s
pTuneos_bin_path=%s
netctl_fold=%s
netchop_path=%s
python ${pTuneos_bin_path}/sm_netMHC_result_parse.py -i ${netmhc_out}/${PREFIX}_all_netmhc.tsv -g ${netmhc_out}/${PREFIX}_all.fasta -o ${netmhc_out} -s ${PREFIX}_all -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
python ${pTuneos_bin_path}/netCTLPAN.py -i ${netmhc_out}/${PREFIX}_all_final_neo_candidate.tsv -d ${driver_gene_path} -o ${netctl_fold} -s ${PREFIX}_all -n ${netchop_path}
'''%(prefix,netmhc_out_fold,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,hla_str,driver_gene_path,pTuneos_bin_path,netctl_out_fold,netchop_path)
	#print str_proc1
	subprocess.call(str_proc1, shell=True, executable='/bin/bash')



def pyclone_annotation(copynumber_fold,somatic_out_fold,prefix,pyclone_fold,netctl_out_fold,pyclone_path,pTuneos_bin_path,logfile_fold,netmhc_out_fold):
	str_proc=r'''
copynumber_fold=%s
somatic_mutation=%s
PREFIX=%s
pyclone=%s
netctl=%s
Pyclone=%s
pTuneos_bin_path=%s
logfile_fold=%s
netmhc_out=%s
python ${pTuneos_bin_path}/sequenza2pyclone.py ${netmhc_out}/${PREFIX}_all_final_neo_candidate.tsv ${somatic_mutation}/${PREFIX}_filter.vcf ${copynumber_fold}/${PREFIX}_segments.txt ${PREFIX} ${pyclone}
TUMOR_CONTENT=`cat ${copynumber_fold}/${PREFIX}_cellularity.txt`
$Pyclone setup_analysis --in_files ${pyclone}/${PREFIX}_sequenza2pyclone.txt --tumour_contents ${TUMOR_CONTENT} --prior major_copy_number --working_dir ${pyclone}
$Pyclone run_analysis --config_file ${pyclone}/config.yaml > ${logfile_fold}/${PREFIX}_pyclone.log 2>&1
$Pyclone build_table --config_file ${pyclone}/config.yaml --out_file ${pyclone}/loci.tsv --table_type loci
python ${pTuneos_bin_path}/neo_pyclone_annotation_vcf.py -n ${netctl}/${PREFIX}_all_netctl_concact.tsv -s ${pyclone}/loci.tsv -o ${netctl} -S ${PREFIX}
'''%(copynumber_fold,somatic_out_fold,prefix,pyclone_fold,netctl_out_fold,pyclone_path,pTuneos_bin_path,logfile_fold,netmhc_out_fold)
	#print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')


def kallisto_expression(raw_fastq_path_first,raw_fastq_path_second,kallisto_path,kallisto_out_fold,prefix,kallisto_cdna_path,logfile_fold,threads,fragment_length,fragment_SD):
	cdna_path_dir = os.path.dirname(kallisto_cdna_path)
	cnd_file_prefix = os.path.splitext(os.path.basename(kallisto_cdna_path))[0]
	kallisto_index_path = cdna_path_dir + '/' + cnd_file_prefix + '.idx'
	if not os.path.exists(kallisto_index_path):
		cmd_kallisto_index = kallisto_path + " index -i " + kallisto_index_path + ' ' + kallisto_cdna_path + ' > ' +  logfile_fold + '/' + prefix + '_kallisto_index.log' + ' 2>&1'
		os.system(cmd_kallisto_index)
	else:
		print "kallisto index already exists. Continue..."
	if os.path.exists(raw_fastq_path_second):
		cmd_kallisto_quant = kallisto_path + " quant -i " + kallisto_index_path + " -t " + str(threads) + " -b 100 -o " + kallisto_out_fold + " " + raw_fastq_path_first + " " + raw_fastq_path_second + ' > ' +  logfile_fold + '/' + prefix + '_kallisto.log' + ' 2>&1'
	else:
		cmd_kallisto_quant = kallisto_path + " quant -i " + kallisto_index_path + " -t " + str(threads) + " -b 100 --single -l " + str(fragment_length) + " -s " + str(fragment_SD) + " -o " + kallisto_out_fold + " " + raw_fastq_path_first +  ' > ' +  logfile_fold + '/' + prefix + '_kallisto.log' + ' 2>&1'
	#print cmd_kallisto_index
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
def InVivoModelAndScore(neo_file,cf_hy_model_9,cf_hy_model_10,cf_hy_model_11,RF_model,neo_model_file,blastp_tmp_file,blastp_out_tmp_file,netMHCpan_pep_tmp_file,netMHCpan_ml_out_tmp_file,iedb_file,blast_db_path):
	iedb_seq=get_iedb_seq(iedb_file)
	hy_xgb_9=joblib.load(cf_hy_model_9)
	hy_xgb_10=joblib.load(cf_hy_model_10)
	hy_xgb_11=joblib.load(cf_hy_model_11)
	#print hy_xgb_9
	data_neo=pd.read_csv(neo_file,header=0,sep='\t')
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
		

		
