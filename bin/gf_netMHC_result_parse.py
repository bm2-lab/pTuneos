#!/usr/bin/python
# -*- coding: UTF-8 -*-
###########gene fusion netMHC result parsing and filter based on binding affinity and FPKM #########
import sys,getopt,os,re
import pandas as pd 


opts,args=getopt.getopt(sys.argv[1:],"hi:t:g:o:b:f:s:l:",["input_netmhc_file","input_genefusion_fasta","ericscript_outfile","out_dir","binding_affinity_cutoff","exp_cutoff","sample_id","hla_str"])
input_netmhc_file=""
input_genefusion_fasta=""
ericscript_outfile =""
out_dir=""
binding_affinity_cutoff=500
exp_cutoff=1
sample_id=""
hla_str=""
USAGE='''usage: python gf_netmhc_result_parse.py -i <input_netmhc_file> -t <input_genefusion_fasta> -g <ericscript_outfile> -o <outdir> -s <sample_id> [option]
		required argument:
			-i | --input_netmhc_file : input file,result from netMHC
			-t | --input_genefusion_fasta : gene fusion derived peptide
			-g | --ericscript_outfile : gene_fusion result from ericscript
			-o | --out_dir : output directory
			-s | --sample_id : sample id
		optional argument:
			-b | --binding_affinity_cutoff : pipetide binding affinity cutoff , default: 500
			-f | --exp_cutoff : expression cutoff, default :10
			-l | --hla_str : hla type string derived from opitype'''
	
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-i","--input_netmhc_file"):
		input_netmhc_file=value
	elif opt in ("-t","--input_genefusion_fasta"):
		input_genefusion_fasta=value
	elif opt in ("-g","--ericscript_outfile"):
		ericscript_outfile =value
	elif opt in ("-o","--out_dir"):
		out_dir =value  
	elif opt in ("-b","--binding_affinity_cutoff"):
		binding_affinity_cutoff=value
	elif opt in ("-f","--exp_cutoff"):
		exp_cutoff=value 
	elif opt in ("-s","sample_id"):
		sample_id=value
	elif opt in ("-l","hla_str"):
		hla_str=value
print binding_affinity_cutoff, exp_cutoff
if (input_netmhc_file =="" or input_genefusion_fasta=="" or ericscript_outfile =="" or out_dir =="" or sample_id=="" or hla_str==""):
	print USAGE
	sys.exit(2)		

######## extract candidate neoantigens####
MT_header = []
for line in open(input_genefusion_fasta,'r'):
    record=line.strip()
    if record.startswith('>'):
        MT_header.append(record[1:])
    else:
        continue

dup_full_header=[]
hla_num=len(hla_str.split(','))
i=0
while i<hla_num:
	for j in range(len(MT_header)):
		dup_full_header.append(MT_header[j])
	i=i+1

with open(input_netmhc_file) as f:
    data = f.read()

nw_data = data.split('-----------------------------------------------------------------------------------\n')

MT_neo = []
for i in range(len(nw_data)):
    if i%4 == 2:
        mt_neo_data = nw_data[i].strip().split('\n')
        MT_neo.append(mt_neo_data)

HLAtype=[]
peptide=[]
fusion_gene=[]
binding_aff=[]
binding_level=[]
binding_level_des=[]
for i in range(len(MT_neo)):
    for j in range(len(MT_neo[i])):
	line=re.split(r'\s+',MT_neo[i][j].strip())
        if line[-1]=="WB" or line[-1]=="SB":
            HLAtype.append(line[1])
            peptide.append(line[2])
            fusion_gene.append(dup_full_header[i])
            binding_aff.append(line[12])
            binding_level.append(line[13])           
            binding_level_des.append(line[-1])
data_eric=pd.read_table(ericscript_outfile,header=0,sep='\t')
fusiongene_eric=[]
for i in range(len(data_eric.GeneName1)):
    fusiongene=data_eric.GeneName1[i]+'>>'+data_eric.GeneName2[i]
    fusiongene_eric.append(fusiongene)
data_eric['Fusion_Gene']=fusiongene_eric
data_netmhc=pd.DataFrame()
data_netmhc['HLA_Gene']=HLAtype
data_netmhc['Peptide']=peptide
data_netmhc['Fusion_Gene']=fusion_gene
data_netmhc['Binding_Aff']=binding_aff
data_netmhc['Binding_Level']=binding_level
data_netmhc['Binding_Level_Des']=binding_level_des
data_merge=pd.merge(data_netmhc,data_eric,on='Fusion_Gene',how='left')
data_select=data_merge.loc[:,['HLA_Gene','Peptide','Fusion_Gene','Binding_Aff','Binding_Level','Binding_Level_Des','GeneExpr_Fused','EricScore']]
data_select.Binding_Aff=data_select.Binding_Aff.astype('float64')
data_select_filter=data_select[(data_select.Binding_Aff<int(binding_affinity_cutoff)) & (data_select.GeneExpr_Fused>int(exp_cutoff))]
data_select_filter_sorted=data_select_filter.sort_values(["Binding_Aff","EricScore"],ascending=[1,0])
data_select_filter_sorted.to_csv(out_dir+'/'+sample_id+"_genefusion_final_neo_candidate.txt",header=1,sep='\t',index=0)

    







            
