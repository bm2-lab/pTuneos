import pandas as pd
import sys,getopt
opts,args=getopt.getopt(sys.argv[1:],"hi:o:s:",["input_seq2HLA_result","out_dir","sample_id"])
input_seq2HLA_result =""
out_dir=""
sample_id=""
USAGE='''
	This script extract hla string from seq2HLA result
	usage: python seq2HLAExt.py -i <input_seq2HLA_result> -o <outdir> -s <sample_id>
		required argument:
			-i | --input_seq2HLA_result : input file result from seq2HLA
			-o | --out_dir : output directory
			-s | --sample_id : sample id
'''
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-i","--input_seq2HLA_result"):
		input_seq2HLA_result=value
	elif opt in ("-o","--out_dir"):
		out_dir =value
	elif opt in ("-s","--sample_id"):
		sample_id =value  
	
#print coverage
if (input_seq2HLA_result =="" or out_dir =="" or sample_id==""):
	print USAGE
	sys.exit(2)
data=pd.read_table(input_seq2HLA_result,header=0)
data_allele=data['Allele 1']
HLA_A='HLA-'+data_allele[0].replace('*','')
HLA_B='HLA-'+data_allele[1].replace('*','')
HLA_C='HLA-'+data_allele[2].replace('*','')
out_str=HLA_A+','+HLA_B+','+HLA_C

f_o=open(out_dir+'/'+sample_id+'_hlatype_result','w')
f_o.write(out_str)
f_o.close()