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
USAGE='''usage: python netCTLPAN.py -i <input_neo_file> -o <outdir> -s <sample_id>
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
if data_neo_fil.empty:
	print "NetMHC filtering result is empty!!"
	sys.exit()
if os.path.exists(out_dir+'/'+sample_id+'_netCLT.txt'):
        os.remove(out_dir+'/'+sample_id+'_netCLT.txt')
for hla,gene,aa,mt_pep in zip(data_neo_fil['HLA_type'],data_neo_fil.Gene,data_neo_fil.AA_change,data_neo_fil.MT_pep):
    line='>'+str(gene) + '\n' + str(mt_pep) + '\n'
    pep_len=len(mt_pep)
    hla_type=hla.replace("*","")
    f=open(out_dir+'/'+sample_id+'_tmp.txt','w')
    f.write(line)
    str_pro='python /home/zhouchi/software/netchop/predict.py --method netctlpan --allele ' + hla_type + ' --length ' +  str(pep_len)+ ' --threshold -99.9 --cleavage_weight 0.225 --tap_weight 0.025 --epitope_threshold 1.0 --noplot ' + out_dir+'/' + sample_id + '_tmp.txt' + ' > '+ out_dir+'/'+sample_id+'_tmp_netCLT.txt'
    f.close()
    subprocess.call(str_pro,shell = True,executable = '/bin/bash')
    flag=1
    while flag:
        if os.path.getsize(out_dir+'/'+sample_id+'_tmp_netCLT.txt')==0:
            subprocess.call(str_pro,shell = True,executable = '/bin/bash')
        else:
            flag=0
    os.system("cat " + out_dir+'/'+sample_id+'_tmp_netCLT.txt' + ">>" + out_dir+'/'+sample_id+'_netCLT.txt')



with open(out_dir+'/'+sample_id+'_netCLT.txt') as f:
    data=f.read()
    
net_res=[line.split('\t') for line in data.strip().split('\n') if line.startswith('1')]
tap_prediction_score=[]
cleavage_prediction_score=[]
combined_prediction_score=[]
for net in net_res:
    tap_prediction_score.append(net[3])
    cleavage_prediction_score.append(net[4])
    combined_prediction_score.append(net[5])
pdata={'tap_prediction_score':tap_prediction_score,'cleavage_prediction_score':cleavage_prediction_score,'combined_prediction_score':combined_prediction_score}

data_pdata=pd.DataFrame(pdata) 
data_con=pd.concat([data_neo_fil,data_pdata],axis=1)
gene_list=data_con.Gene.drop_duplicates()
f_drivergene=open('/home/zhouchi/database/Annotation/DriveGene.tsv','r')
drivergene_list=[]
for line in f_drivergene:
	drivergene_list.append(line.strip())
whether_driver_gene=[]
for i in range(len(data_con.Gene)):
	if data_con.Gene[i] in drivergene_list:
		whether_driver_gene.append('TRUE')
	else:
		whether_driver_gene.append('FALSE')
data_con['DriverGene_Lable'] = whether_driver_gene
#data_con_drop=data_con.drop_duplicates(subset=['#Position','HLA_type','Gene','Mutation','MT_pep','WT_pep'])
data_con_drop=data_con.drop_duplicates(subset=['#Position','HLA_type','Gene','MT_pep','WT_pep'])
data_con_sort=data_con_drop.sort_values(['MT_Binding_Aff','fold_change'],ascending=[1, 1])
data_has_change_aa=data_con_sort[data_con_sort.MT_pep!=data_con_sort.WT_pep]
data_has_change_aa=data_has_change_aa.reset_index()
del data_has_change_aa['index']
data_has_change_aa["contain_X"]=['X' not in data_has_change_aa.MT_pep[i] for i in range(len(data_has_change_aa.MT_pep))]
data_out=data_has_change_aa[data_has_change_aa.contain_X==True]
data_out.fillna(0.5)
data_out.to_csv(out_dir+'/'+sample_id+'_netctl_concact.tsv',sep='\t',header=1,index=0)	 
if os.path.exists(out_dir+'/'+sample_id+'_tmp.txt'):    
	os.remove(out_dir+'/'+sample_id+'_tmp.txt')
if os.path.exists(out_dir+'/'+sample_id+'_netCLT.txt'):
	os.remove(out_dir+'/'+sample_id+'_netCLT.txt')
if os.path.exists(out_dir+'/'+sample_id+'_tmp_netCLT.txt'):
	os.remove(out_dir+'/'+sample_id+'_tmp_netCLT.txt')