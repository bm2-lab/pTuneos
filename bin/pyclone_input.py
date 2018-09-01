#!/usr/bin/python
# -*- coding: UTF-8 -*-
###########process varscan result into pyclone input#########
import pandas as pd
import sys,getopt

opts,args=getopt.getopt(sys.argv[1:],"hn:i:s:c:o:S:",["neoantigen_file","vep_input_file","snv_mutect2_file","copynumber_file","out_dir","SampleID"])
neoantigen_file=""
vep_input_file=""
snv_mutect2_file=""
copynumber_file=""
output_dir=""
sample_id=""

USAGE='''usage: python pyclone_input.py -n <neoantigen_file> -i <vep_input_file> -s <snv_mutect2_file> -c <copynumber_file> -o <outdir> -S <sample_id> [option]"
		required argument:
			-n | --neoantigen_file: 
			-i | --vep_input_file : snv vep input file
			-s | --snv_mutect2_file : snp file result from mutect2
			-c | --copynumber_file : copynumber file result from varscan
			-o | --out_dir : output_directory
			-S | --SampleID : sample id

'''
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-n","--neoantigen_file"):
		neoantigen_file=value
	elif opt in ("-i","--vep_input_file"):
		vep_input_file=value
	elif opt in ("-s","--snv_mutect2_file"):
		snv_mutect2_file=value
	elif opt in ("-c","--copynumber_file"):
		copynumber_file=value
	elif opt in ("-o","--out_dir"):
		output_dir=value  
	elif opt in ("-S","--SampleID"):
		sample_id=value

if (neoantigen_file=="" or vep_input_file=="" or snv_mutect2_file=="" or copynumber_file=="" or output_dir=="" or sample_id==""):
	print USAGE
	sys.exit(2)
gene_dic={}	 ###gene--postion dictionary
vep_pos_list=[]#####postion information in VEP
vep_gene_list=[]####gene information in VEP
f_vep=open(vep_input_file)
for ele in f_vep:
	if ele.startswith('#'):
		continue
	else:
		vep_pos=ele.strip().split('\t')[1]
		extra=ele.strip().split('\t')[-1]
		#vep_gene=extra.split(';')[2].split('=')[1]
		element=extra.split(';')
		for ele in element:
			sub_ele=ele.split('=')
			if sub_ele[0]=="SYMBOL":
				vep_gene=sub_ele[1]
	if vep_gene not in vep_gene_list:
		vep_pos_list.append(vep_pos)
		vep_gene_list.append(vep_gene)

f_vep.close()
#print len(vep_pos_list),len(vep_gene_list)
for i in range(len(vep_pos_list)):
	gene_dic[vep_gene_list[i]]=vep_pos_list[i]
#print gene_dic,len(gene_dic)
neo_chr_pos=[]
gene_name=[]
for ele in open(neoantigen_file,'r'):
	if ele.startswith('#'):
		continue
	else:
		gene=ele.strip().split('\t')[2]
		n_c_r=gene_dic[gene]
	if gene not in gene_name:
		gene_name.append(gene)
		neo_chr_pos.append(n_c_r)

print len(neo_chr_pos)
#test=[]
#test1=[]
#for ele in neo_chr_pos:
#	if ele not in test1:
#		test1.append(ele)
#		continue
#	else:
#		test.append(ele)
#print test
#print len(vep_pos_list),len(vep_gene_list)
#print vep_pos_list
#print vep_gene_list
#chr_pos_list=[]
#gene_list=[]
#for i in range(len(vep_pos_list)):
#	chr_pos=vep_pos_list[i].split(':')[0]+':'+vep_pos_list[i].split(':')[1]
#	chr_pos_list.append(chr_pos)
#	gene_list.append(vep_gene_list[i])

chr_list=[]
pos_list=[]
ref_allele=[]
var_allele=[]
tumor_reads1=[]
tumor_reads2=[]
snp_chr_pos=[]
f_snp=open(snv_mutect2_file,'r')
for ele in f_snp:
	if ele.startswith('#'):
		continue
	else:
		line=ele.strip().split('\t')
		chro_pos=line[0]+':'+line[1]
		c_l=line[0]
		p_l=line[1]
		r_a=line[2]
		a_a=line[3]
		format_information=line[8].split(':')
		alt_read_1_index=format_information.index('ALT_F1R2')
		alt_read_2_index=format_information.index('ALT_F2R1')
		ref_read_1_index=format_information.index('REF_F1R2')
		ref_read_2_index=format_information.index('REF_F2R1')
		tumor_info_list=line[9].split(':')
		print tumor_info_list
		t_r1=int(tumor_info_list[ref_read_1_index])+int(tumor_info_list[ref_read_2_index])
		t_r2=int(tumor_info_list[alt_read_1_index])+int(tumor_info_list[alt_read_2_index])
		#t_r1=line[8]
		#t_r2=line[9]
	if chro_pos in neo_chr_pos:
		snp_chr_pos.append(chro_pos)
		chr_list.append(c_l)
		pos_list.append(p_l)
		ref_allele.append(r_a)
		var_allele.append(a_a)
		tumor_reads1.append(t_r1)
		tumor_reads2.append(t_r2)
####finnaly, the list neo_chr_pos is duplicated,however,this will not affect our analysis##
f_snp.close()
data_snp=pd.DataFrame()
data_snp['chrom']=chr_list
data_snp['position']=pos_list
data_snp['ref']=ref_allele
data_snp['var']=var_allele
data_snp['tumor_reads1']=tumor_reads1
data_snp['tumor_reads2']=tumor_reads2
#print data_snp
data_copynumber=pd.read_table(copynumber_file,header=0,sep='\t')
data_allele_count=data_snp[['chrom','position','ref','var','tumor_reads1','tumor_reads2']]
data_cn_count=data_copynumber[['chromosome','start.pos','end.pos','CNt','A','B']]
range_dic={}
for i in range(len(data_cn_count.chromosome)):
	if data_cn_count.chromosome[i] not in range_dic.keys():
		range_dic[data_cn_count.chromosome[i]]=[]
		range_dic[data_cn_count.chromosome[i]].append([data_cn_count['start.pos'][i],data_cn_count['end.pos'][i],data_cn_count.CNt[i],data_cn_count.A[i],data_cn_count.B[i]])
	else:
		range_dic[data_cn_count.chromosome[i]].append([data_cn_count['start.pos'][i],data_cn_count['end.pos'][i],data_cn_count.CNt[i],data_cn_count.A[i],data_cn_count.B[i]])
CNT_list=[]
CNA_list=[]
CNB_list=[]
#print range_dic
for i in range(len(data_allele_count.chrom)):
	chr_str=data_allele_count.chrom[i]
	pos=data_allele_count.position[i]
	flag=0
	for ele in range_dic[chr_str]:
		#flag=0
		start_pos=ele[0]
		end_pos=ele[1]
		CNT=ele[2]
		CNA=ele[3]
		CNB=ele[4]
		if long(pos)>=start_pos and long(pos)<=end_pos:
			CNT_list.append(CNT)
			CNA_list.append(CNA)
			CNB_list.append(CNB)
			flag=1
			break
	if flag==0:
		CNT_list.append(2)
		CNA_list.append(2)
		CNB_list.append(0)

			
###Genotype####
#print CNA_list,len(CNA_list)
#print CNB_list,len(CNB_list)
mutation_case_name=sample_id
data_allele_count['variant_case']=mutation_case_name
data_allele_count['normal_cn']=2
data_allele_count['total_copy_number']=CNT_list
data_allele_count['major_cn']=CNA_list
#print CNA_list,len(CNA_list)
data_allele_count['minor_cn']=CNB_list
#print CNB_list,len(CNB_list)
Genotype=[]
#j=0
for i in range(len(data_allele_count.major_cn)):
	if data_allele_count.major_cn[i]>0 and data_allele_count.minor_cn[i]>0:
		gt="AB"
	elif data_allele_count.major_cn[i]>0 and data_allele_count.minor_cn[i]==0:
		gt="BB"
	elif data_allele_count.major_cn[i]==0 and data_allele_count.minor_cn[i]>0:
		gt="AA"
	else:
		print "Ops,this cause by chrosome Y total copynumber undetecting!"
#		print data_allele_count.major_cn[i],data_allele_count.minor_cn[i]
#		j=j+1

		
	Genotype.append(gt)
#print j
data_allele_count['genotype']=Genotype
#####mutaton ID####
MutationID=[]
for i in range(len(data_allele_count.chrom)):
	mt_id=':'.join([mutation_case_name, data_allele_count.genotype[i], data_allele_count.chrom[i], str(data_allele_count.position[i])])
	MutationID.append(mt_id)
data_allele_count['mutation_id']=MutationID
VariantFreq=[]
for i in range(len(data_allele_count.chrom)):
	vf=float(data_allele_count.tumor_reads2[i])/(float(data_allele_count.tumor_reads1[i])+float(data_allele_count.tumor_reads2[i]))
	VariantFreq.append(vf)
data_allele_count['variant_freq']=VariantFreq
data_pyclone_input=data_allele_count[['mutation_id','tumor_reads1','tumor_reads2','normal_cn','minor_cn','major_cn','variant_case','variant_freq','genotype']]
data_pyclone_input=data_pyclone_input.rename(columns={'tumor_reads1':'ref_counts','tumor_reads2':'var_counts'})
data_pyclone_input.to_csv(output_dir+'/'+sample_id+'_pyclone_input.tsv',sep='\t',index=0)



