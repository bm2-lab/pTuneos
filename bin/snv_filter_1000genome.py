import pandas as pd
import numpy as np
import sys,getopt,os
import re
#####prepare fasta format input file for netMHC######
opts,args=getopt.getopt(sys.argv[1:],"hi:g:o:s:",["input_vcf_file","input_snv_1000G_file","out_dir","sample_id"])
input_vcf_file=""
input_snv_1000G_file =""
out_dir=""
sample_id=""
USAGE='''
	This script convert VEP result to fasta format file for netMHC
	usage: python animo_acid_prepare.py -i <input_vcf_file> -g <input_snv_1000G_file> -o <outdir> -s <sample_id>
		required argument:
			-i | --input_vcf_file : input vcf file
			-g | --1000G_pos : input 1000G snp position file 
			-o | --out_dir : output directory
			-s | --sample_id : sample id
'''
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-i","--input_vcf_file"):
		input_vcf_file=value
	elif opt in ("-g","--input_snv_1000G_file"):
		input_snv_1000G_file=value
	elif opt in ("-o","--out_dir"):
		out_dir =value
	elif opt in ("-s","--sample_id"):
		sample_id =value  
	
#print coverage
if (input_vcf_file =="" or out_dir =="" or sample_id==""):
	print USAGE
	sys.exit(2)	
snv_1000G=[]
for line in open(input_snv_1000G_file):
	snv_1000G.append(line.strip())
set_snv=set(snv_1000G)
f_out=open(out_dir+'/'+sample_id+"_SNVs_filter_1000.vcf",'w')

for line in open(input_vcf_file):
	if line.startswith('#'):
		f_out.write(line)
		continue
	else:
		record=line.strip().split('\t')
		pos=record[0]+':'+record[1]
		if pos not in set_snv:
			f_out.write(line)
		else:
			continue
f_out.close()






