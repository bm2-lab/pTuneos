#!/usr/bin/python
# -*- coding: UTF-8 -*-
###########netMHC result parsing and filter based on binding affinity and FPKM #########
import sys,getopt,os,subprocess
import pandas as pd 
import time,os


opts,args=getopt.getopt(sys.argv[1:],"hi:d:v:n:o:s:",["input_vcf_file","tumor_depth_cutoff","tumor_vaf_cutoff","normal_vaf_cutoff","out_dir","sample_id"])
input_vcf_file =""
tumor_depth_cutoff=7
tumor_vaf_cutoff=0.05
normal_vaf_cutoff=0.03
out_dir=""
sample_id=""
USAGE='''usage: python snv_filter.py -i <input_neo_file> -d <tumor_depth_cutoff> -v <tumor_vaf_cutoff> -n <normal_vaf_cutoff> -o <outdir> -s <sample_id>
		required argument:
			-i | --input_vcf_file : input vcf file
			-d | --tumor_depth_cutoff : tumor depth cutoff
			-v | --tumor_vaf_cutoff :tumor vaf cutoff
			-n | --normal_vaf_cutoff : normal vaf cutoff
			-o | --out_dir : output directory
			-s | --sample_id : sample id
'''
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-i","--input_vcf_file"):
		input_vcf_file=value
	elif opt in ("-d","--tumor_depth_cutoff"):
		tumor_depth_cutoff=value
	elif opt in ("-v","--tumor_vaf_cutoff"):
		tumor_vaf_cutoff=value
	elif opt in ("-n","--normal_vaf_cutoff"):
		normal_vaf_cutoff=value
	elif opt in ("-o","--out_dir"):
		out_dir =value  
	elif opt in ("-s","sample_id"):
		sample_id=value
#print coverage
if (input_vcf_file =="" or out_dir =="" or sample_id==""):
	print USAGE
	sys.exit(2)


f_filter=open(out_dir+'/'+sample_id+"_filter.vcf",'w')

for line in open(input_vcf_file):
	if line.startswith("#"):
		f_filter.write(line)
	else:
		record=line.strip().split('\t')
		tumor_info=record[9]
		normal_info=record[10]
		tumor_depth=tumor_info.split(':')[1].split(',')[1]
		tumor_vaf=tumor_info.split(':')[2]
		normal_vaf=normal_info.split(':')[2]
		if int(tumor_depth)>=float(tumor_depth_cutoff) and float(tumor_vaf)>float(tumor_vaf_cutoff) and float(normal_vaf)<float(normal_vaf_cutoff):
			f_filter.write(line)

f_filter.close()