import sys,getopt
opts,args=getopt.getopt(sys.argv[1:],"hi:o:s:",["indel_file","out_dir","sample_id"])
indel_file ="" 
out_dir=""
sample_id=""
USAGE='''
	This script convert deletion VCF derived VEP result to fasta format file for netMHC
	usage: python deletion2fasta.py -i <input_vep_file> -o <outdir> -s <sample_id>
		required argument:
			-i | --indel_file : input file,result from VEP
			-o | --out_dir : output directory
			-s | --sample_id : sample id
'''
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-i","--indel_file"):
		indel_file=value
	elif opt in ("-o","--out_dir"):
		out_dir =value
	elif opt in ("-s","--sample_id"):
		sample_id =value  
	
#print coverage
if (indel_file =="" or out_dir =="" or sample_id==""):
	print USAGE
	sys.exit(2)	
chrom_list=[]
start_list=[]
end_list=[]
ref_alt_list=[]

for ele in open(indel_file):
	if ele.startswith('chrom'):
		continue
	else:
		line=ele.strip().split('\t')
		#print line
		if line[12]=='Somatic':
			chrom_list.append(line[0])
			if line[3].startswith('+'):
				start=int(line[1])+1
				end=int(line[1])
				ref_alt='-'+'/'+line[3][1:]
			else:
				start=int(line[1])+1
				end=int(line[1])+len(line[3][1:])
				ref_alt=line[3][1:]+'/'+'-'
			start_list.append(start)
			end_list.append(end)
			ref_alt_list.append(ref_alt)
		
#print len(chrom_list),len(ref_alt_list)
f_o=open(out_dir+'/'+sample_id+'_varscan_indel.vcf','w')
for i in range(len(chrom_list)):
	f_o.write(chrom_list[i]+'\t'+str(start_list[i])+'\t'+str(end_list[i])+'\t'+ref_alt_list[i]+'\n')
f_o.close()
