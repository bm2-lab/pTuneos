#!/usr/bin/python
import sys,getopt
opts,args=getopt.getopt(sys.argv[1:],"hi:o:s:",["input_optitype_result","out_dir","sample_id"])
input_optitype_result =""
out_dir=""
sample_id=""
USAGE='''
	This script extract hla string from optitype result
	usage: python optitype_ext.py -i <input_optitype_result> -o <outdir> -s <sample_id>
		required argument:
			-i | --input_optitype_result : input file result from optitype
			-o | --out_dir : output directory
			-s | --sample_id : sample id
'''
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-i","--input_optitype_result"):
		input_optitype_result=value
	elif opt in ("-o","--out_dir"):
		out_dir =value
	elif opt in ("-s","--sample_id"):
		sample_id =value  
	
#print coverage
if (input_optitype_result =="" or out_dir =="" or sample_id==""):
	print USAGE
	sys.exit(2)

with open(input_optitype_result) as f:
	data=f.read()
line_list=data.strip().split('\n')[1].split('\t')[1:-2]
line_list_modify=[]
for line in line_list:
	line_modify='HLA-'+line.replace('*','')
	line_list_modify.append(line_modify)

line_out=','.join(line_list_modify)
f_o=open(out_dir+'/'+sample_id+'_optitype_hla_type','w')
f_o.write(line_out)
f_o.close()