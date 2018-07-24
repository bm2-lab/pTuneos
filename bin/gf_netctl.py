#!/usr/bin/python
# -*- coding: UTF-8 -*-
###########netMHC result parsing and filter based on binding affinity and FPKM #########
import sys,getopt,os,subprocess
import pandas as pd 
import time,os


opts,args=getopt.getopt(sys.argv[1:],"hi:o:s:",["input_neo_file","out_dir","sample_id"])
input_neo_file =""
out_dir=""
sample_id=""
USAGE='''usage: python gf_netctl.py -i <input_neo_file> -o <outdir> -s <sample_id>
		required argument:
			-i | --input_neo_file : input file,result from netMHC parse
			-o | --out_dir : output directory
			-s | --sample_id : sample id
'''
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-i","--input_neo_file"):
		input_neo_file=value
	elif opt in ("-o","--out_dir"):
		out_dir =value  
	elif opt in ("-s","sample_id"):
		sample_id=value
#print coverage
if (input_neo_file =="" or out_dir =="" or sample_id==""):
	print USAGE
	sys.exit(2)


data_neo_fil = pd.read_table(input_neo_file,header=0,sep='\t')
if os.path.exists(out_dir+'/'+sample_id+'_netCLT.txt'):
    	os.remove(out_dir+'/'+sample_id+'_netCLT.txt')
for hla,gene,mt_pep in zip(data_neo_fil['HLA_Gene'],data_neo_fil["Fusion_Gene"],data_neo_fil["Peptide"]):
    fusion_gene=gene.replace(">>","_")
    line='>'+str(fusion_gene) + '\n' + str(mt_pep) + '\n'
    pep_len=len(mt_pep)
    #hla_type=hla[0:7]+':'+hla[7:]
    hla_type=hla.replace("*","")
    f=open(out_dir+'/'+sample_id+'_tmp.txt','w')
    f.write(line)
    str_pro='python /home/zhouchi/software/netchop/predict.py --method netctlpan --allele ' + hla_type + ' --length ' +  str(pep_len)+ ' --threshold -99.9 --cleavage_weight 0.225 --tap_weight 0.025 --epitope_threshold 1.0 --noplot ' + out_dir+'/' + sample_id + '_tmp.txt' + ' >> '+ out_dir+'/'+sample_id+'_netCLT.txt'
    print str_pro
    f.close()
    subprocess.call(str_pro,shell = True,executable = '/bin/bash')



with open(out_dir+'/'+sample_id+'_netCLT.txt') as f:
    data=f.read()
    
net_res=[line.split('\t') for line in data.strip().split('\n') if line.startswith('1')]
tap_prediction_score=[]
cleavage_prediction_score=[]
combined_prediction_score=[]
for net in net_res:
    tap_prediction_score.append(net[2])
    cleavage_prediction_score.append(net[3])
    combined_prediction_score.append(net[4])
pdata={'tap_prediction_score':tap_prediction_score,'cleavage_prediction_score':cleavage_prediction_score,'combined_prediction_score':combined_prediction_score}

data_pdata=pd.DataFrame(pdata) 
data_con=pd.concat([data_neo_fil,data_pdata],axis=1)
gene_list=data_con["Fusion_Gene"].drop_duplicates()
f_drivergene=open('/home/zhouchi/database/Annotation/GoldDriverGene','r')
drivergene_list=[]
for line in f_drivergene:
	drivergene_list.append(line)
whether_driver_gene=[]
for i in range(len(data_con["Fusion_Gene"])):
	if data_con["Fusion_Gene"][i] in drivergene_list:
		whether_driver_gene.append('TRUE')
	else:
		whether_driver_gene.append('FALSE')
data_con['DriverGene_Lable'] = whether_driver_gene
data_con.to_csv(out_dir+'/'+sample_id+'_netctl_concact.txt',sep='\t',header=1,index=0)	 
if os.path.exists(out_dir+'/'+sample_id+'_tmp.txt'):    
	os.remove(out_dir+'/'+sample_id+'_tmp.txt')
if os.path.exists(out_dir+'/'+sample_id+'_netCLT.txt'):
	os.remove(out_dir+'/'+sample_id+'_netCLT.txt')