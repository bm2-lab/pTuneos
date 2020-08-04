#!/usr/bin/python
# -*- coding: UTF-8 -*-
###########process varscan result into pyclone input#########
import pandas as pd
import numpy as np
import sys,getopt

opts,args=getopt.getopt(sys.argv[1:],"hn:s:o:S:",["neoantigen_file","pyclone_loci","out_dir","SampleID"])
neoantigen_file=""
pyclone_loci=""
output_dir=""
sample_id=""

USAGE='''usage: python pyclone_input.py -n <neoantigen_file> -s <pyclone_loci> -o <outdir> -S <sample_id> [option]"
		required argument:
			-n | --neoantigen_file: 
			-s | --pyclone_loci : pyclone_loci file result from pyclone
			-o | --out_dir : output_directory
			-S | --SampleID : sample id
'''
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-n","--neoantigen_file"):
		neoantigen_file=value
	elif opt in ("-s","--pyclone_loci"):
		pyclone_loci=value
	elif opt in ("-o","--out_dir"):
		out_dir=value  
	elif opt in ("-S","--SampleID"):
		sample_id=value


#print coverage
if (neoantigen_file=="" or pyclone_loci=="" or out_dir=="" or sample_id==""):
	print USAGE
	sys.exit(2)

pos_cp_dic={}
pos_vaf_dic={}
data_loci=pd.read_csv(pyclone_loci,header=0,sep='\t')
for i in range(len(data_loci.mutation_id)):
	pos=data_loci.mutation_id[i].split(':')[1]+":"+data_loci.mutation_id[i].split(':')[2]
	pos_cp_dic[pos]=data_loci.cellular_prevalence[i]
	pos_vaf_dic[pos]=data_loci.variant_allele_frequency[i]

VAF_list=[]
CP_list=[]
data_neo=pd.read_csv(neoantigen_file,header=0,sep='\t')
for i in range(len(data_neo["#Position"])):
	neo_pos=data_neo["#Position"][i].split(":")[1].split("-")
	if len(neo_pos)==1:
		CP=pos_cp_dic[data_neo["#Position"][i]]
		VAF=pos_vaf_dic[data_neo["#Position"][i]]
	else:
		neo_pos_new=data_neo["#Position"][i].split(":")[0]+":"+str(long(neo_pos[0])-1)
		CP=pos_cp_dic[neo_pos_new]
		VAF=pos_vaf_dic[neo_pos_new]
	VAF_list.append(VAF)
	CP_list.append(CP)

data_neo["cellular_prevalence"]=CP_list
data_neo["variant_allele_frequency"]=VAF_list

del data_neo["cleavage_prediction_score"]
del data_neo["tap_prediction_score"]
#print data_neo
data_fill_na=data_neo.fillna(0)
####HLA:0201 and sort by aff_rank
data_drop=data_fill_na.drop_duplicates(subset=["Gene","MT_pep","WT_pep"])
data_sort=data_drop.sort_values(["MT_Binding_Aff"],ascending=True)
del data_sort["Transcript_name"]
#print data_sort
data_sort.to_csv(out_dir+'/'+sample_id+'_pyclone_neo.tsv',header=1,sep='\t',index=0)


















