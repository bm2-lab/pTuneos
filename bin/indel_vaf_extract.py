#!/usr/bin/python
# -*- coding: UTF-8 -*-
###########get indel vaf information from varscan and strelka#########
import pandas as pd
import numpy as np
import sys,getopt

opts,args=getopt.getopt(sys.argv[1:],"hn:i:s:o:S:",["neoantigen_file","varscan_indel_file","strelka_indel_file","out_dir","SampleID"])
neoantigen_file=""
varscan_indel_file=""
strelka_indel_file=""
output_dir=""
sample_id=""
USAGE='''usage: python pyclone_input.py -n <neoantigen_file> -i <varscan_indel_file> -s <strelka_indel_file> -o <outdir> -S <sample_id> [option]"
		required argument:
			-n | --neoantigen_file: 
			-i | --varscan_indel : varscan indel input file
			-s | --strelka_indel : strelka indel file
			-o | --out_dir : output_directory
			-S | --SampleID : sample id
'''
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-n","--neoantigen_file"):
		neoantigen_file=value
	elif opt in ("-i","--varscan_indel"):
		varscan_indel_file=value
	elif opt in ("-s","--strelka_indel"):
		strelka_indel_file=value
	elif opt in ("-o","--out_dir"):
		out_dir=value  
	elif opt in ("-S","--SampleID"):
		sample_id=value

#print coverage
if (neoantigen_file=="" or varscan_indel_file=="" or strelka_indel_file=="" or out_dir=="" or sample_id==""):
	print USAGE
	sys.exit(2)



indel_vaf_dic={}
for line in open(varscan_indel_file):
	if not line.startswith("chrom"):
		record=line.strip().split('\t')
		if record[12]=="Somatic":
			pos=record[0]+":"+record[1]
			tumor_vaf=record[10]
			indel_vaf_dic[pos]=tumor_vaf


for line in open(strelka_indel_file):
	if not line.startswith("#"):
		record=line.strip().split('\t')
		if record[6]=="PASS":
			pos=record[0]+":"+record[1]
			tumor_info=record[10]
			print record
			if float(tumor_info.split(':')[2].split(',')[0]) != 0:
				tumor_vaf=float(tumor_info.split(':')[3].split(',')[0])/float(tumor_info.split(':')[2].split(',')[0])
			else:
				tumor_vaf=0.10
			indel_vaf_dic[pos]=tumor_vaf


f_out=open(out_dir+'/'+sample_id+"_neo_indel_vaf.txt",'w')
#indel_vaf=[]
for line in open(neoantigen_file):
	if line.startswith('#'):
		f_out.write(line.strip()+"\t"+"variant_allele_frequency"+"\n")
	else:
		record=line.strip().split('\t')
		pos=int(record[0].split('-')[0].split(':')[1])-1
		chr_pos=record[0].split('-')[0].split(':')[0]+":"+str(pos)
		for key in indel_vaf_dic.keys():
			key_pos=key.split(":")[1]
			if abs(int(key_pos)-pos)<=2:
				#indel_vaf.append(indel_vaf_dic[chr_pos])
				if str(indel_vaf_dic[key])[-1]=="%":
					out_indel_vaf=float(indel_vaf_dic[key][:-1])/100
				else:
					out_indel_vaf=indel_vaf_dic[key]
				f_out.write(line.strip()+'\t'+str(out_indel_vaf)+'\n')


f_out.close()






