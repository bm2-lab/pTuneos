import subprocess
import os,sys,time
import os,sys,time
import multiprocessing
import shutil 
import subprocess
import pandas as pd
import math
from pyper import *
import numpy as np
import math
from sklearn import preprocessing
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from sklearn.semi_supervised import label_propagation
import itertools
from scipy import linalg
import matplotlib as mpl
from sklearn import mixture
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from imblearn.under_sampling import RandomUnderSampler
from imblearn.metrics import classification_report_imbalanced
from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from math import log, exp
import pandas as pd
import math	
import numpy as np
from scipy import interp
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
import xgboost as xgb
from xgboost.sklearn import XGBClassifier
from sklearn import cross_validation, metrics
from sklearn.grid_search import GridSearchCV
import matplotlib.pylab as plt
from sklearn.model_selection import train_test_split
from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from math import log, exp
import subprocess
from sklearn.ensemble import RandomForestClassifier  
from sklearn.preprocessing import StandardScaler
from imblearn.over_sampling import SMOTE
from collections import Counter
from sklearn.cross_validation import cross_val_score
a=26
k=4.86936
M=1. #default concentration of mutant peptides
W=1. #default concentration of wildtype peptides

WEPS=0.0003
HYDROPHOBIC_RESIDUES="AILMFWYV"
WEIRD_RESIDUES="CGP"
hydro_score={"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,"K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,"T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3}
def get_iedb_seq(iedb_file):
	iedb_seq=[]
	for line in open(iedb_file):
		if line.startswith(">"):
			continue
		else:
			iedb_seq.append(line.strip())
	return iedb_seq
def VCF_process(prefix,vcf_file,somatic_out_fold,vcftools_path,vep_path,vep_cache_path,netmhc_out_path,tumor_depth_cutoff,tumor_vaf_cutoff,normal_vaf_cutoff):
	cmd_mutation_filter='grep ' + "\'^#\|chr[1-9]\{0,1\}[0-9XY]\\{0,1\\}\\b\'" + ' ' + vcf_file + ' > ' + somatic_out_fold + '/' + prefix + '_' + 'mutect2_filter.vcf'
	print cmd_mutation_filter
	os.system(cmd_mutation_filter)
	cmd_vcftools_snv=vcftools_path + " --vcf " + vcf_file + " --remove-filtered-all --remove-indels --recode --recode-INFO-all --out " + somatic_out_fold + '/' + prefix + '_'+ 'SNVs_only'
	cmd_vcftools_indel=vcftools_path + " --vcf " + vcf_file + " --remove-filtered-all --keep-only-indels --recode --recode-INFO-all --out " + somatic_out_fold + '/' + prefix + '_'+ 'INDELs_only'
	print cmd_vcftools_snv
	print cmd_vcftools_indel
	os.system(cmd_vcftools_snv)
	os.system(cmd_vcftools_indel)
	cmd_snv_filter="python ${iTuNES_BIN_PATH}/snv_filter.py -i " + somatic_out_fold + '/' + prefix + '_'+ 'SNVs_only.recode.vcf' + " -d " + str(tumor_depth_cutoff) + " -v " + str(tumor_vaf_cutoff) + " -n " + str(normal_vaf_cutoff) + " -o " + somatic_out_fold + " -s " + prefix
	print cmd_snv_filter
	os.system(cmd_snv_filter)
	cmd_vep=vep_path + " -i " + somatic_out_fold + '/' + prefix + '_'+ 'SNVs_filter.vcf' + " --cache --dir " + vep_cache_path + " --dir_cache " + vep_cache_path + " --force_overwrite --canonical --symbol -o STDOUT --offline | filter_vep --ontology --filter \"CANONICAL is YES and Consequence is missense_variant\" -o " + somatic_out_fold + '/' + prefix + '_'+ 'snv_vep_ann.txt' + " --force_overwrite"
	print cmd_vep
	os.system(cmd_vep)
	cmd_vep_snv_all=vep_path + " -i " + somatic_out_fold + '/' + prefix + '_'+ 'SNVs_filter.vcf' + " --cache --dir " + vep_cache_path + " --dir_cache " + vep_cache_path + " --force_overwrite --canonical --symbol -o STDOUT --offline | filter_vep --ontology --filter \"Consequence is missense_variant\" -o " + somatic_out_fold + '/' + prefix + '_'+ 'snv_vep_ann_all.txt' + " --force_overwrite"
	print cmd_vep_snv_all
	os.system(cmd_vep_snv_all)
	cmd_vep_indel=vep_path + " -i " + somatic_out_fold + '/' + prefix + '_'+ 'INDELs_only.recode.vcf' + " --cache --dir " + vep_cache_path + " --dir_cache " + vep_cache_path + " --force_overwrite --canonical --symbol -o STDOUT --offline | filter_vep --ontology --filter \"Consequence is missense_variant\" -o " + somatic_out_fold + '/' + prefix + '_'+ 'mutect_indel_vep_ann.txt' + " --force_overwrite"
	print cmd_vep_indel
	os.system(cmd_vep_indel)
	cmd_snv="python ${iTuNES_BIN_PATH}/snv2fasta.py -i " + somatic_out_fold + '/' + prefix + '_'+ 'snv_vep_ann.txt' + " -o " + netmhc_out_path + " -s " + prefix
	print cmd_snv
	os.system(cmd_snv)

def netMHCpan(fasta_file,hla_str,netmhc_out_file,out_dir,split_num,netMHCpan_path,tmp_dir,peptide_length):
	str_proc=r'''
set -x
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
echo ${arr1[@]}
OLD_IFS="$IFS" 
IFS=","
arr2=(${hla_str})
IFS="$OLD_IFS" 
for s in ${arr2[@]}
do
{
	echo $s
	for file_l in ${arr1[@]}
	do
	{
		echo ${file_l}
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
set +x
'''%(fasta_file,hla_str,netmhc_out_file,out_dir,split_num,netMHCpan_path,tmp_dir,peptide_length)
	subprocess.call(str_proc, shell=True, executable='/bin/bash')



def snv_neo(snv_fasta_file,hla_str,snv_netmhc_out_file,netmhc_out_fold,split_num,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_out_fold,netMHCpan_path,peptide_length):
	netMHCpan(snv_fasta_file,hla_str,snv_netmhc_out_file,netmhc_out_fold,split_num,netMHCpan_path,'tmp_snv',peptide_length)
	str_proc1=r'''
PREFIX=%s
netmhc_out=%s
Exp_file=%s
Binding_Aff_Fc_Cutoff=%d
Binding_Aff_Cutoff=%d
Fpkm_Cutoff=%d
hla_str=%s
netctl_fold=%s
python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i ${netmhc_out}/${PREFIX}_snv_netmhc.txt -g ${netmhc_out}/${PREFIX}_snv.fasta -o ${netmhc_out} -s ${PREFIX}_snv -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
python ${iTuNES_BIN_PATH}/netCTLPAN.py -i ${netmhc_out}/${PREFIX}_snv_final_neo_candidate.txt -o ${netctl_fold} -s ${PREFIX}_snv
'''%(prefix,netmhc_out_fold,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,hla_str,netctl_out_fold)
	print str_proc1
	subprocess.call(str_proc1, shell=True, executable='/bin/bash')


def indel_neo(indel_fasta_file,somatic_out_fold,hla_str,netmhc_out_file,split_num,netMHCpan_path,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_out_fold,netmhc_out_fold,peptide_length):
	count = 0
	for index, line in enumerate(open(somatic_out_fold + "/" + prefix + "_mutect_indel_vep_ann.txt",'r')):
		count += 1
	if count == 45:
		print "No indel sites were detected!"
	else:
		str_proc1=r'''
	PREFIX=%s
	somatic_fold=%s
	netmhc_out=%s
	python ${iTuNES_BIN_PATH}/varscandel2fasta.py -i ${somatic_fold}/${PREFIX}_indel_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}
	python ${iTuNES_BIN_PATH}/varscanins2fasta.py -i ${somatic_fold}/${PREFIX}_indel_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}
	cat ${netmhc_out}/${PREFIX}_del.fasta > ${netmhc_out}/${PREFIX}_indel.fasta
	cat ${netmhc_out}/${PREFIX}_ins.fasta >> ${netmhc_out}/${PREFIX}_indel.fasta
	'''%(prefix,somatic_out_fold,netmhc_out_fold)
		#subprocess.call(str_proc1, shell=True, executable='/bin/bash')
		netMHCpan(indel_fasta_file,hla_str,netmhc_out_file,netmhc_out_fold,split_num,netMHCpan_path,"tmp_indel",peptide_length)
		str_proc2=r'''
	PREFIX=%s
	netmhc_out=%s
	Exp_file=%s
	Binding_Aff_Fc_Cutoff=%d
	Binding_Aff_Cutoff=%d
	Fpkm_Cutoff=%d
	hla_str=%s
	netctl_fold=%s
	python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i ${netmhc_out}/${PREFIX}_indel_netmhc.txt -g ${netmhc_out}/${PREFIX}_indel.fasta -o ${netmhc_out} -s ${PREFIX} -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff}
	python ${iTuNES_BIN_PATH}/netCTLPAN.py -i ${netmhc_out}/${PREFIX}_indel_final_neo_candidate.txt -o ${netctl_fold} -s ${PREFIX}
	'''%(prefix,netmhc_out_fold,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,hla_str,netctl_out_fold)	
		#subprocess.call(str_proc2, shell=True, executable='/bin/bash')


def pyclone_annotation(somatic_out_fold,copynumber_profile,tumor_cellularity,prefix,pyclone_fold,netctl_out_fold,coverage,pyclone_path):
	str_proc=r'''
somatic_mutation=%s
copynumber_profile=%s
TUMOR_CONTENT=%f
PREFIX=%s
pyclone=%s
netctl=%s
COVERAGE=%d
Pyclone=%s
python ${iTuNES_BIN_PATH}/pyclone_input.py -n ${netctl}/${PREFIX}_snv_netctl_concact.txt -i ${somatic_mutation}/${PREFIX}_snv_vep_ann_all.txt -s ${somatic_mutation}/${PREFIX}_SNVs_only.recode.vcf -c ${copynumber_profile} -o ${pyclone} -S ${PREFIX} -C ${COVERAGE}
$Pyclone setup_analysis --in_files ${pyclone}/${PREFIX}_pyclone_input.tsv --tumour_contents $TUMOR_CONTENT --prior major_copy_number --working_dir ${pyclone}
$Pyclone run_analysis --config_file ${pyclone}/config.yaml
$Pyclone build_table --config_file ${pyclone}/config.yaml --out_file ${pyclone}/loci.tsv --table_type loci
python ${iTuNES_BIN_PATH}/neo_pyclone_annotation.py -n ${netctl}/${PREFIX}_snv_netctl_concact.txt -i ${somatic_mutation}/${PREFIX}_snv_vep_ann_all.txt -s ${pyclone}/loci.tsv -o ${netctl} -S ${PREFIX}
'''%(somatic_out_fold,copynumber_profile,tumor_cellularity,prefix,pyclone_fold,netctl_out_fold,coverage,pyclone_path)
	print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')

def hydro_vector(pep):
	hydro_score={"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,"K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,"T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3}
	hydrophobicity_vector=[]
	pep_list=list(pep)
	pep_len=len(pep_list)
	for pep in pep_list:
		hydrophobicity_vector.append(hydro_score[pep.upper()])
	return hydrophobicity_vector
def Train_hydrophobicity_xgboost(postive_pep_file,negative_pep_file):
	pos_pep_list=[]
	for line in open(postive_pep_file):
		if line.startswith(">"):
			continue
		else:
			record=line.strip()
			if len(record)==9:
				pos_pep_list.append(record)
	neg_pep_list=[]
	for line in open(negative_pep_file):
		if line.startswith(">"):
			continue
		else:
			record=line.strip()
			if len(record)==9:
				neg_pep_list.append(record)
	pos_pep_hydro_list=[hydro_vector(p) for p in pos_pep_list]#[0:len(neg_pep_all)]
	neg_pep_hydro_list=[hydro_vector(p) for p in neg_pep_list]
	X=pos_pep_hydro_list+neg_pep_hydro_list
	X_scale = StandardScaler().fit_transform(X)
	y=[1]*len(pos_pep_hydro_list)+[0]*len(neg_pep_hydro_list)
	nm1 = NearMiss(random_state=0, version=1)
	X_resampled_nm1, y_resampled_nm1 = nm1.fit_sample(X, y)
	clf_mlp = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(3), random_state=0)
	clf_mlp.fit(X_resampled_nm1, y_resampled_nm1)
	#hydro_score=clf_mlp.predict_prob([hydro_vector[pep]])
	return clf_mlp

#def Train_bioactivity_rf():
	


def aligner(seq1,seq2):
	matrix = matlist.blosum62
	gap_open = -11
	gap_extend = -1
	aln = pairwise2.align.localds(seq1.upper(), seq2.upper(), matrix,gap_open,gap_extend)
	return aln

def logSum(v):
	ma=max(v)
	return log(sum(map(lambda x: exp(x-ma),v)))+ma



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



def GetNmerPositivePep(n,mhc_pos_file):
	pep_list=[]
	for line in open(mhc_pos_file):
		if line.startswith(">"):
			continue
		else:
			record=line.strip()
			if len(record)==n:
				pep_list.append(record)
	return pep_list


def GetNmerNegativePep(n,mhc_neg_file):
	pep_list=[]
	for line in open(mhc_neg_file):
		if line.startswith(">"):
			continue
		else:
			record=line.strip()
			if len(record)==n:
				pep_list.append(record)
	return pep_list

def hydro_vector(pep):
	hydrophobicity_vector=[]
	pep_list=list(pep)
	pep_len=len(pep_list)
	for pep in pep_list:
		hydrophobicity_vector.append(hydro_score[pep.upper()])

	return hydrophobicity_vector

def getXY(n,mhc_pos_file,mhc_neg_file):
	pos_pep=GetNmerPositivePep(n,mhc_pos_file)
	neg_pep=GetNmerNegativePep(n,mhc_neg_file)
	####get hydrophobicity list
	pos_pep_hydro_list=[hydro_vector(p) for p in pos_pep]#[0:len(neg_pep_all)]
	neg_pep_hydro_list=[hydro_vector(p) for p in neg_pep]
	#print len(pos_pep_hydro_list)
	#print len(neg_pep_hydro_list)
	######transform into array
	#pos_pep_hydro_array=np.array(pos_pep_hydro_list)
	#neg_pep_hydro_array=np.array(neg_pep_hydro_list)
	X=pos_pep_hydro_list+neg_pep_hydro_list
	y=[1]*len(pos_pep_hydro_list)+[0]*len(neg_pep_hydro_list)
	X_array=np.asarray(X)
	y_array=np.asarray(y)
	return X_array,y_array
def get_hydro_model(mhc_pos_file,mhc_neg_file):
	xgb_9 = XGBClassifier(
	 learning_rate =0.1,
	 n_estimators=208,
	 max_depth=5,
	 min_child_weight=1,
	 gamma=0,
	 subsample=0.8,
	 colsample_bytree=0.8,
	 objective= 'binary:logistic',
	 n_jobs=4,
	 scale_pos_weight=1,
	 random_state=27)
	#[207]  train-auc:0.95222+0.000337282   test-auc:0.852987+0.00777271
	#AUC Score (Train): 0.941597
	#AUC Score (Test): 0.768922

	xgb_10 = XGBClassifier(
	 learning_rate =0.1,
	 n_estimators=165,
	 max_depth=16,
	 min_child_weight=1,
	 gamma=0,
	 subsample=0.8,
	 colsample_bytree=0.8,
	 objective= 'binary:logistic',
	 n_jobs=4,
	 scale_pos_weight=1,
	 random_state=27)
	#[164]  train-auc:1+0   test-auc:0.900267+0.00955733
	#AUC Score (Train): 1.000000
	#AUC Score (Test): 0.802474
	xgb_11 = XGBClassifier(
	 learning_rate =0.1,
	 n_estimators=123,
	 max_depth=5,
	 min_child_weight=1,
	 gamma=0.1,
	 subsample=0.8,
	 colsample_bytree=0.8,
	 objective= 'binary:logistic',
	 n_jobs=4,
	 scale_pos_weight=1,
	 random_state=27)
	X_array_9,y_array_9=getXY(9,mhc_pos_file,mhc_neg_file)
	hy_xgb_9=xgb_9.fit(X_array_9,y_array_9)
	X_array_10,y_array_10=getXY(10,mhc_pos_file,mhc_neg_file)
	hy_xgb_10=xgb_10.fit(X_array_10,y_array_10)
	X_array_11,y_array_11=getXY(11,mhc_pos_file,mhc_neg_file)
	hy_xgb_11=xgb_11.fit(X_array_11,y_array_11)
	return hy_xgb_9,hy_xgb_10,hy_xgb_11
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
	print str_blastp_pro
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
	print human_homolog_pep
	f=open(netMHCpan_pep_tmp_file,'w')
	f.write(human_homolog_pep+'\n')
	f.close()
	str_netMHCpan_ml_pro='netMHCpan -p ' + netMHCpan_pep_tmp_file + ' -a ' + hla_type_in + ' > ' + netMHCpan_ml_out_tmp_file
	print str_netMHCpan_ml_pro
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
	print str_netMHCpan_ml_pro
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
def InVivoModelAndScoreSNV(mhc_pos_file,mhc_neg_file,neo_file,model_train_file,neo_model_file,blastp_tmp_file,blastp_out_tmp_file,netMHCpan_pep_tmp_file,netMHCpan_ml_out_tmp_file,iedb_file,blast_db_path,immunogenicity_gmm_all_score_ranking,immunogenicity_gmm_pos_score_ranking,gmm_classification_file,immunogenicity_bioactive_score_ranking):
	iedb_seq=get_iedb_seq(iedb_file)
	hy_xgb_9,hy_xgb_10,hy_xgb_11=get_hydro_model(mhc_pos_file,mhc_neg_file)
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
			print "Oh no!!"
			print line
			print len(line)
			hydrophobicity_score.append(0.5)
			R=calculate_R(line,iedb_seq)	
			Recognition_score.append(R)
	paired_similarity_score=[]
	homolog_similaity_score=[]
	#####paired similarity and homolog similarity########
	for M_P,N_P,H_P in zip(data_neo.MT_pep,data_neo.WT_pep,Homolog_pep):
		print M_P,N_P,H_P
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
	data_train = pd.read_table(model_train_file,header=0,sep='\t')
	data_train_dropna=data_train.dropna()
	target = 'response'
	HPcol = 'hydrophobicity_score'
	Rcol = 'Recognition_score'
	SScol = 'self_sequence_similarity'
	ELcol = 'EL_dissimilarity'
	predictors = [x for x in data_train_dropna.columns if x not in [target,'Pep_len',ELcol,'EL_ht_rank']]
	X_train=data_train_dropna[predictors].values
	X_train_scale = StandardScaler().fit_transform(X_train)
	y_train=data_train_dropna[target].values
	print 'Original dataset shape {}'.format(Counter(y_train))
	sm=SMOTE(k_neighbors=4,kind='borderline1',random_state=42)
	X_res, y_res = sm.fit_sample(X_train,y_train)
	print 'Resampled dataset shape {}'.format(Counter(y_res))
	rf0 = RandomForestClassifier(oob_score=True, random_state=10)  
	print rf0
	rf0.fit(X_res, y_res)
	dneo_predprob = rf0.predict_proba(df_neo.values)[:,1]
	print dneo_predprob
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
	data_neo_out_sort.to_csv(neo_model_file,sep='\t',header=1,index=0)
	del data_neo_out_sort["contain_X"]
def InVivoModelAndScoreINDEL(mhc_pos_file,mhc_neg_file,neo_file,model_train_file,neo_model_file,blastp_tmp_file,blastp_out_tmp_file,netMHCpan_pep_tmp_file,netMHCpan_ml_out_tmp_file,iedb_file,blast_db_path,immunogenicity_gmm_all_score_ranking,immunogenicity_gmm_pos_score_ranking,gmm_classification_file,immunogenicity_bioactive_score_ranking):
	iedb_seq=get_iedb_seq(iedb_file)
	hy_xgb_9,hy_xgb_10,hy_xgb_11=get_hydro_model(mhc_pos_file,mhc_neg_file)
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
			print "Oh no!!"
			print line
			print len(line)
			hydrophobicity_score.append(0.5)
			R=calculate_R(line,iedb_seq)	
			Recognition_score.append(R)
	paired_similarity_score=[]
	homolog_similaity_score=[]
	#####paired similarity and homolog similarity########
	for M_P,N_P,H_P in zip(data_neo.MT_pep,data_neo.WT_pep,Homolog_pep):
		print M_P,N_P,H_P
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
	data_train = pd.read_table(model_train_file,header=0,sep='\t')
	data_train_dropna=data_train.dropna()
	target = 'response'
	HPcol = 'hydrophobicity_score'
	Rcol = 'Recognition_score'
	SScol = 'self_sequence_similarity'
	ELcol = 'EL_dissimilarity'
	predictors = [x for x in data_train_dropna.columns if x not in [target,'Pep_len',ELcol,'EL_ht_rank']]
	X_train=data_train_dropna[predictors].values
	X_train_scale = StandardScaler().fit_transform(X_train)
	y_train=data_train_dropna[target].values
	print 'Original dataset shape {}'.format(Counter(y_train))
	sm=SMOTE(k_neighbors=4,kind='borderline1',random_state=42)
	X_res, y_res = sm.fit_sample(X_train,y_train)
	print 'Resampled dataset shape {}'.format(Counter(y_res))
	rf0 = RandomForestClassifier(oob_score=True, random_state=10)  
	print rf0
	rf0.fit(X_res, y_res)
	dneo_predprob = rf0.predict_proba(df_neo.values)[:,1]
	print dneo_predprob
	data_neo["model_pro"]=dneo_predprob
	neoantigen_infor = []
	for i in range(len(data_neo.Gene)):
		neo_infor = data_neo["Gene"][i]+'_'+data_neo["AA_change"][i]+'_'+data_neo["MT_pep"][i]+'_'+data_neo["WT_pep"][i]
		neoantigen_infor.append(neo_infor)
	data_neo["neoantigen_infor"] = neoantigen_infor
	data_neo_out_sort=data_neo.sort_values(['model_pro'],ascending=[0])
	data_neo_out_sort.to_csv(neo_model_file,sep='\t',header=1,index=0)
	f_EL_rank_wt=lambda x:1-(1/(1+math.pow(math.e,5*(float(x)-2))))/2
	f_EL_rank_mt=lambda x:1/(1+math.pow(math.e,5*(float(x)-2)))
	EL_mt_rank_score=data_neo.MT_Binding_EL.apply(f_EL_rank_mt)
	EL_wt_rank_score=data_neo.WT_Binding_EL.apply(f_EL_rank_wt)
	bioactive_score=[data_neo.Hydrophobicity_score[i]*data_neo.Recognition_score[i]*data_neo.Self_sequence_similarity[i]*EL_mt_rank_score[i]*EL_wt_rank_score[i] for i in range(len(data_neo.MT_Binding_EL))]
	data_neo["bioactive_score"]=bioactive_score
	data_neo_sortedby_bioactive=data_neo.sort_values(["bioactive_score"],ascending=False)
	data_neo_sortedby_bioactive.to_csv(immunogenicity_bioactive_score_ranking,header=1,index=0,sep='\t')
	k=1
	f_TPM=lambda x:math.tanh(x/k)
	allele_frequency_score=data_neo.variant_allele_frequency
	netchop_score=data_neo.combined_prediction_score
	try:
		tpm_score=data_neo.tpm.apply(f_TPM)
		immuno_effect_score=[tpm_score[i]*allele_frequency_score[i]*netchop_score[i]*data_neo.Hydrophobicity_score[i]*data_neo.Recognition_score[i]*data_neo.Self_sequence_similarity[i]*EL_mt_rank_score[i]*EL_wt_rank_score[i] for i in range(len(data_neo.MT_Binding_EL))]
		data_feature_select=pd.DataFrame()
		data_feature_select["neoantigen_infor"]=neoantigen_infor
		data_feature_select["EL_mt_rank_score"]=EL_mt_rank_score
		data_feature_select["EL_wt_rank_score"]=EL_wt_rank_score
		data_feature_select["tpm_score"]=tpm_score
		#data_feature_select["tpm_normal_score"]=tpm_normal_score
		data_feature_select["allele_frequency_score"]=allele_frequency_score
		data_feature_select["netchop_score"]=netchop_score
		data_feature_select["hydrophobicity_score"]=data_neo.Hydrophobicity_score
		data_feature_select["Recogonition"]=data_neo.Recognition_score
		data_feature_select["Self_sequence_similarity"]=data_neo.Self_sequence_similarity
		data_feature_select["immuno_effect_score"]=immuno_effect_score
		X_neo = data_feature_select.values[:,1:9]
		####pca on non-standardized data of all 6 feature
		pca = PCA(n_components=2).fit(X_neo)
		X_neo_pca = pca.transform(X_neo)
		#print "PCA result on non-standardized data"
		#print "explained_variance_ratio",pca.explained_variance_ratio_
		#print "explained_variance",pca.explained_variance_
		#print pca.n_components_
		#print pca.components_
		###plot GMM ellipse
		# Fit a Gaussian mixture with EM using five components
		gmm = mixture.GaussianMixture(n_components=2, covariance_type='full',n_init=5).fit(X_neo_pca)
		plot_results(X_neo_pca, gmm.predict(X_neo_pca), gmm.means_, gmm.covariances_, 0,'Gaussian Mixture',gmm_classification_file)
		#plot_results(X_neo_pca, gmm.predict(X_neo_pca), gmm.means_, gmm.covars_, 0, 'Gaussian Mixture')
		predict_label=gmm.predict(X_neo_pca)
		predict_prob=gmm.predict_proba(X_neo_pca)
		predict_positive_prob=[]
		for i in range(len(predict_prob)):
			pos_prob=predict_prob[i][1]
			predict_positive_prob.append(pos_prob)
		data_feature_select["gmm_label"]=predict_label
		data_gmm_filter_1=data_feature_select[data_feature_select["gmm_label"]==1]
		data_gmm_filter_1_scoreAve=data_gmm_filter_1.immuno_effect_score.sum()/len(data_gmm_filter_1.immuno_effect_score)
		data_gmm_filter_0=data_feature_select[data_feature_select["gmm_label"]==0]
		data_gmm_filter_0_scoreAve=data_gmm_filter_0.immuno_effect_score.sum()/len(data_gmm_filter_0.immuno_effect_score)
		if data_gmm_filter_0_scoreAve > data_gmm_filter_1_scoreAve:
			data_gmm_filter_0.gmm_label=1
			data_gmm_filter_1.gmm_label=0
		else:
			pass
		data_labeled=pd.concat([data_gmm_filter_1,data_gmm_filter_0])
		data_label_sorted=data_labeled.sort_values(["immuno_effect_score"],ascending=False)
		data_label_sorted.to_csv(immunogenicity_gmm_all_score_ranking,header=1,index=0,sep='\t')
		data_gmm_positive=data_label_sorted[data_label_sorted["gmm_label"]==1]
		data_gmm_positive.to_csv(immunogenicity_gmm_pos_score_ranking,header=1,index=0,sep='\t')
	except AttributeError,e:
		print e
		immuno_effect_score=[allele_frequency_score[i]*netchop_score[i]*data_neo.Hydrophobicity_score[i]*data_neo.Recognition_score[i]*data_neo.Self_sequence_similarity[i]*EL_mt_rank_score[i]*EL_wt_rank_score[i] for i in range(len(data_neo.MT_Binding_EL))]
		data_feature_select=pd.DataFrame()
		data_feature_select["neoantigen_infor"]=neoantigen_infor
		data_feature_select["EL_mt_rank_score"]=EL_mt_rank_score
		data_feature_select["EL_wt_rank_score"]=EL_wt_rank_score
		data_feature_select["allele_frequency_score"]=allele_frequency_score
		data_feature_select["netchop_score"]=netchop_score
		data_feature_select["hydrophobicity_score"]=data_neo.Hydrophobicity_score
		data_feature_select["Recogonition"]=data_neo.Recognition_score
		data_feature_select["Self_sequence_similarity"]=data_neo.Self_sequence_similarity
		data_feature_select["immuno_effect_score"]=immuno_effect_score
		X_neo = data_feature_select.values[:,1:8]
		####pca on non-standardized data of all 6 feature
		pca = PCA(n_components=2).fit(X_neo)
		X_neo_pca = pca.transform(X_neo)
		#print "PCA result on non-standardized data"
		#print "explained_variance_ratio",pca.explained_variance_ratio_
		#print "explained_variance",pca.explained_variance_
		#print pca.n_components_
		#print pca.components_
		###plot GMM ellipse
		# Fit a Gaussian mixture with EM using five components
		gmm = mixture.GaussianMixture(n_components=2, covariance_type='full',n_init=5).fit(X_neo_pca)
		plot_results(X_neo_pca, gmm.predict(X_neo_pca), gmm.means_, gmm.covariances_, 0,'Gaussian Mixture',gmm_classification_file)
		#plot_results(X_neo_pca, gmm.predict(X_neo_pca), gmm.means_, gmm.covars_, 0, 'Gaussian Mixture')
		predict_label=gmm.predict(X_neo_pca)
		predict_prob=gmm.predict_proba(X_neo_pca)
		predict_positive_prob=[]
		for i in range(len(predict_prob)):
			pos_prob=predict_prob[i][1]
			predict_positive_prob.append(pos_prob)
		data_feature_select["gmm_label"]=predict_label
		data_gmm_filter_1=data_feature_select[data_feature_select["gmm_label"]==1]
		data_gmm_filter_1_scoreAve=data_gmm_filter_1.immuno_effect_score.sum()/len(data_gmm_filter_1.immuno_effect_score)
		data_gmm_filter_0=data_feature_select[data_feature_select["gmm_label"]==0]
		data_gmm_filter_0_scoreAve=data_gmm_filter_0.immuno_effect_score.sum()/len(data_gmm_filter_0.immuno_effect_score)
		if data_gmm_filter_0_scoreAve > data_gmm_filter_1_scoreAve:
			data_gmm_filter_0.gmm_label=1
			data_gmm_filter_1.gmm_label=0
		else:
			pass
		data_labeled=pd.concat([data_gmm_filter_1,data_gmm_filter_0])
		data_label_sorted=data_labeled.sort_values(["immuno_effect_score"],ascending=False)
		data_label_sorted.to_csv(immunogenicity_gmm_all_score_ranking,header=1,index=0,sep='\t')
		data_gmm_positive=data_label_sorted[data_label_sorted["gmm_label"]==1]
		data_gmm_positive.to_csv(immunogenicity_gmm_pos_score_ranking,header=1,index=0,sep='\t')











