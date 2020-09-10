###process somatic mutation
import sys
import pandas as pd

neo_file=sys.argv[1]
mutation_file=sys.argv[2]
segment_copynumber_file=sys.argv[3]
sample_name=sys.argv[4]
out_dir=sys.argv[5]


data_neo=pd.read_csv(neo_file,header=0,sep='\t')

data_neo_dup=data_neo.drop_duplicates(["#Position"])

chr_pos_list_raw=list(data_neo_dup["#Position"])
chr_pos_list=[]
#print len(chr_pos_list)
for line in chr_pos_list_raw:
	pos_list_sub=line.split(":")[1].split("-")
	if len(pos_list_sub)==1:
		new_line=line
	else:
		new_line=line.split(":")[0]+":"+str(long(pos_list_sub[0])-1)
	chr_pos_list.append(new_line)




mutation_list=[]
chr_list=[]
pos_list=[]
ref_allele=[]
var_allele=[]
ref_reads=[]
alt_reads=[]
snp_chr_pos=[]
VAF=[]
f_mutation=open(mutation_file,'r')
for ele in f_mutation:
	if ele.startswith('#'):
		continue
	else:
		line=ele.strip().split('\t')
		chro_pos_0=line[0]+':'+line[1]
		chro_pos_1=line[0]+':'+str(long(line[1])-1)
		chro_pos_2=line[0]+':'+str(long(line[1])+1)
		if line[0]!="chrMT" and (chro_pos_0 in chr_pos_list):
			chro_pos=line[0]+':'+line[1]
			chr_name=line[0]
			pos_loc=line[1]
			ref_n=line[2]
			alt_n=line[3]
			mut_id=sample_name+":"+chro_pos
			tumor_read_info=line[9].split(':')[1].split(',')
			alt_count=int(tumor_read_info[1])
			ref_count=int(tumor_read_info[0])
			vaf=float(line[9].split(':')[2])
			chr_list.append(chr_name)
			pos_list.append(pos_loc)
			ref_allele.append(ref_n)
			var_allele.append(alt_n)
			ref_reads.append(ref_count)
			alt_reads.append(alt_count)
			snp_chr_pos.append(chro_pos)
			mutation_list.append(mut_id)
			VAF.append(vaf)
		elif line[0]!="chrMT" and (chro_pos_1 in chr_pos_list):
			chro_pos=line[0]+':'+str(long(line[1])-1)
			chr_name=line[0]
			pos_loc=line[1]
			ref_n=line[2]
			alt_n=line[3]
			mut_id=sample_name+":"+chro_pos
			tumor_read_info=line[9].split(':')[1].split(',')
			alt_count=int(tumor_read_info[1])
			ref_count=int(tumor_read_info[0])
			vaf=float(line[9].split(':')[2])
			chr_list.append(chr_name)
			pos_list.append(pos_loc)
			ref_allele.append(ref_n)
			var_allele.append(alt_n)
			ref_reads.append(ref_count)
			alt_reads.append(alt_count)
			snp_chr_pos.append(chro_pos)
			mutation_list.append(mut_id)
			VAF.append(vaf)
		elif line[0]!="chrMT" and (chro_pos_2 in chr_pos_list):
			chro_pos=line[0]+':'+str(long(line[1])+1)
			chr_name=line[0]
			pos_loc=line[1]
			ref_n=line[2]
			alt_n=line[3]
			mut_id=sample_name+":"+chro_pos
			tumor_read_info=line[9].split(':')[1].split(',')
			alt_count=int(tumor_read_info[1])
			ref_count=int(tumor_read_info[0])
			vaf=float(line[9].split(':')[2])
			chr_list.append(chr_name)
			pos_list.append(pos_loc)
			ref_allele.append(ref_n)
			var_allele.append(alt_n)
			ref_reads.append(ref_count)
			alt_reads.append(alt_count)
			snp_chr_pos.append(chro_pos)
			mutation_list.append(mut_id)
			VAF.append(vaf)
		else:
			continue

import pandas as pd 
data_snp=pd.DataFrame()
data_snp["mutation_id"]=mutation_list
data_snp['chrom']=chr_list
data_snp['position']=pos_list
data_snp['ref_counts']=ref_reads
data_snp['var_counts']=alt_reads



data_copynumber=pd.read_csv(segment_copynumber_file,header=0,sep='\t')
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
for i in range(len(data_snp.chrom)):
	chr_str=data_snp.chrom[i]
	pos=data_snp.position[i]
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

data_snp['normal_cn']=2
data_snp['minor_cn']=CNB_list
data_snp['major_cn']=CNA_list
data_snp['variant_case']="test"
data_snp["variant_freq"]=VAF

genotype_list=[]

for i in range(len(data_snp.chrom)):
	if int(data_snp.minor_cn[i])==int(data_snp.major_cn[i]):
		genotype="AB"
	else:
		genotype="BB"
	genotype_list.append(genotype)
data_snp["genotype"]=genotype_list

del data_snp["chrom"]
del data_snp["position"]

data_snp.to_csv(out_dir+"/"+sample_name+"_sequenza2pyclone.txt",header=1,sep="\t",index=0)












