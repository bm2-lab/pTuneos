import os,sys,time
import multiprocessing
import shutil 
import subprocess
import pandas as pd
import math
import numpy as np
from sklearn import preprocessing
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from sklearn.semi_supervised import label_propagation
import matplotlib as mpl
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
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
hydro_score={"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,"K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,"T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3}
def get_iedb_seq(iedb_file):
	iedb_seq=[]
	for line in open(iedb_file):
		if line.startswith(">"):
			continue
		else:
			iedb_seq.append(line.strip())
	return iedb_seq



def VCF_process(prefix,vcf_file,somatic_out_fold,vcftools_path,vep_path,vep_cache_path,netmhc_out_path,tumor_depth_cutoff,tumor_vaf_cutoff,normal_vaf_cutoff,pTuneos_bin_path,human_peptide_path,logfile_fold):
	cmd_mutation_filter='grep ' + "\'^#\|chr[1-9]\{0,1\}[0-9XY]\\{0,1\\}\\b\'" + ' ' + vcf_file + ' > ' + somatic_out_fold + '/' + prefix + '_' + 'mutect2_filter.vcf'
	#print cmd_mutation_filter
	os.system(cmd_mutation_filter)
	cmd_vcftools_snv=vcftools_path + " --vcf " + vcf_file + " --remove-filtered-all --remove-indels --recode --recode-INFO-all --out " + somatic_out_fold + '/' + prefix + '_'+ 'SNVs_only > ' + logfile_fold + '/' + prefix + '_vcftools_snv.log 2>&1'
	cmd_vcftools_indel=vcftools_path + " --vcf " + vcf_file + " --remove-filtered-all --keep-only-indels --recode --recode-INFO-all --out " + somatic_out_fold + '/' + prefix + '_'+ 'INDELs_only > ' + logfile_fold + '/' + prefix + '_vcftools_indel.log 2>&1'
	#print cmd_vcftools_snv
	#print cmd_vcftools_indel
	os.system(cmd_vcftools_snv)
	os.system(cmd_vcftools_indel)
	cmd_snv_filter="python " + pTuneos_bin_path + "/snv_filter.py -i " + somatic_out_fold + '/' + prefix + '_'+ 'SNVs_only.recode.vcf' + " -d " + str(tumor_depth_cutoff) + " -v " + str(tumor_vaf_cutoff) + " -n " + str(normal_vaf_cutoff) + " -o " + somatic_out_fold + " -s " + prefix
	#print cmd_snv_filter
	os.system(cmd_snv_filter)
	cmd_vep=vep_path + " -i " + somatic_out_fold + '/' + prefix + '_'+ 'SNVs_filter.vcf' + " --cache --dir " + vep_cache_path + " --dir_cache " + vep_cache_path + " --force_overwrite --canonical --symbol -o STDOUT --offline | filter_vep --ontology --filter \"CANONICAL is YES and Consequence is missense_variant\" -o " + somatic_out_fold + '/' + prefix + '_'+ 'snv_vep_ann.txt' + " --force_overwrite"
	#print cmd_vep
	os.system(cmd_vep)
	cmd_vep_snv_all=vep_path + " -i " + somatic_out_fold + '/' + prefix + '_'+ 'SNVs_filter.vcf' + " --cache --dir " + vep_cache_path + " --dir_cache " + vep_cache_path + " --force_overwrite --canonical --symbol -o STDOUT --offline | filter_vep --ontology --filter \"Consequence is missense_variant\" -o " + somatic_out_fold + '/' + prefix + '_'+ 'snv_vep_ann_all.txt' + " --force_overwrite"
	#print cmd_vep_snv_all
	os.system(cmd_vep_snv_all)
	cmd_vep_indel=vep_path + " -i " + somatic_out_fold + '/' + prefix + '_'+ 'INDELs_only.recode.vcf' + " --cache --dir " + vep_cache_path + " --dir_cache " + vep_cache_path + " --force_overwrite --canonical --symbol -o STDOUT --offline | filter_vep --ontology --filter \"CANONICAL is YES and Consequence is coding_sequence_variant\" -o " + somatic_out_fold + '/' + prefix + '_'+ 'mutect_indel_vep_ann.txt' + " --force_overwrite"
	#print cmd_vep_indel
	os.system(cmd_vep_indel)
	cmd_vep_indel_all=vep_path + " -i " + somatic_out_fold + '/' + prefix + '_'+ 'INDELs_only.recode.vcf' + " --cache --dir " + vep_cache_path + " --dir_cache " + vep_cache_path + " --force_overwrite --canonical --symbol -o STDOUT --offline | filter_vep --ontology --filter \"Consequence is coding_sequence_variant\" -o " + somatic_out_fold + '/' + prefix + '_'+ 'mutect_indel_vep_ann_all.txt' + " --force_overwrite"
	#print cmd_vep_indel_all
	os.system(cmd_vep_indel_all)
	cmd_snv="python " + pTuneos_bin_path + "/snv2fasta.py -i " + somatic_out_fold + '/' + prefix + '_'+ 'snv_vep_ann.txt' + " -o " + netmhc_out_path + " -s " + prefix + " -p " + human_peptide_path
	#print cmd_snv
	os.system(cmd_snv)


def VCF_process_new(prefix,vcf_file,somatic_out_fold,vcftools_path,vep_path,vep_cache_path,netmhc_out_path,tumor_depth_cutoff,tumor_vaf_cutoff,normal_vaf_cutoff,pTuneos_bin_path,human_peptide_path,logfile_fold,REFERENCE):
	cmd_mutation_filter='grep ' + "\'^#\|chr[1-9]\{0,1\}[0-9XY]\\{0,1\\}\\b\'" + ' ' + vcf_file + ' > ' + somatic_out_fold + '/' + prefix + '_' + 'mutect2_filter.vcf'
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




def pyclone_annotation(somatic_out_fold,copynumber_profile,tumor_cellularity,prefix,pyclone_fold,netctl_out_fold,pyclone_path,pTuneos_bin_path,logfile_fold,netmhc_out_fold):
	str_proc=r'''
somatic_mutation=%s
copynumber_profile=%s
TUMOR_CONTENT=%f
PREFIX=%s
pyclone=%s
pTuneos_bin_path=%s
netctl=%s
Pyclone=%s
logfile_fold=%s
netmhc_out=%s
python ${pTuneos_bin_path}/sequenza2pyclone.py ${netmhc_out}/${PREFIX}_all_final_neo_candidate.tsv ${somatic_mutation}/${PREFIX}_filter.vcf ${copynumber_profile} ${PREFIX} ${pyclone}
#python ${pTuneos_bin_path}/pyclone_input.py -n ${netctl}/${PREFIX}_snv_netctl_concact.tsv -i ${somatic_mutation}/${PREFIX}_snv_vep_ann_all.txt -s ${somatic_mutation}/${PREFIX}_SNVs_only.recode.vcf -c ${copynumber_profile} -o ${pyclone} -S ${PREFIX}
$Pyclone setup_analysis --in_files ${pyclone}/${PREFIX}_sequenza2pyclone.txt --tumour_contents $TUMOR_CONTENT --prior major_copy_number --working_dir ${pyclone}
$Pyclone run_analysis --config_file ${pyclone}/config.yaml > ${logfile_fold}/${PREFIX}_pyclone.log 2>&1
$Pyclone build_table --config_file ${pyclone}/config.yaml --out_file ${pyclone}/loci.tsv --table_type loci
python ${pTuneos_bin_path}/neo_pyclone_annotation_vcf.py -n ${netctl}/${PREFIX}_all_netctl_concact.tsv -s ${pyclone}/loci.tsv -o ${netctl} -S ${PREFIX}
#python ${pTuneos_bin_path}/neo_pyclone_annotation.py -n ${netctl}/${PREFIX}_snv_netctl_concact.tsv -i ${somatic_mutation}/${PREFIX}_snv_vep_ann_all.txt -s ${pyclone}/loci.tsv -o ${netctl} -S ${PREFIX}
'''%(somatic_out_fold,copynumber_profile,tumor_cellularity,prefix,pyclone_fold,pTuneos_bin_path,netctl_out_fold,pyclone_path,logfile_fold,netmhc_out_fold)
	#print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')

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
			print "Oh no!!"
			#print line
			#print len(line)
			hydrophobicity_score.append(0.5)
			R=calculate_R(line,iedb_seq)	
			Recognition_score.append(R)
	paired_similarity_score=[]
	homolog_similaity_score=[]
	#####paired similarity and homolog similarity########
	for M_P,N_P,H_P in zip(data_neo.MT_pep,data_neo.WT_pep,Homolog_pep):
		#print M_P,N_P,H_P
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
	#print dneo_predprob
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












