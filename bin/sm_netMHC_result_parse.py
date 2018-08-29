#!/usr/bin/python
# -*- coding: UTF-8 -*-
###########netMHC result parsing and filter based on binding affinity and FPKM #########
import sys,getopt,os
import pandas as pd 


opts,args=getopt.getopt(sys.argv[1:],"hi:g:e:o:b:a:f:s:l:",["input_netmhc_file","input_fasta","expression_fpkm_file","out_dir","binding_affinity_cutoff","binding_affinity_foldchange_cutoff","fpkm_cutoff","sample_id","hla_str"])
input_netmhc_file=""
input_fasta=""
expression_fpkm_file=""
out_dir=""
binding_affinity_cutoff=500
binding_affinity_foldchange_cutoff=1
fpkm_cutoff=1
sample_id=""
hla_str=""
USAGE='''usage: python netMHC_result_parse.py -i <input_netmhc_file> -g <input_fasta> -o <outdir> -s <sample_id> [option]
		required argument:
			-i | --input_netmhc_file : input file,result from netMHC
			-g | --input_fasta : input fasta file for netMHC
			-o | --out_dir : output directory
			-s | --sample_id : sample id
		optional argument:
			-e | --expression_fpkm_file : expression profile with fpkm value
			-b | --binding_affinity_cutoff : pipetide binding affinity cutoff , default: 500
			-a | --binding_affinity_foldchange_cutoff : pipetide binding affinity fold change between mutate peptide and wild type peptide cutoff , default: 1
			-f | --fpkm_cutoff : FPKM value cutoff, default :1
			-l | --hla_str : hla type string derived from opitype'''
	
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-i","--input_netmhc_file"):
		input_netmhc_file=value
	elif opt in ("-g","--input_fasta"):
		input_fasta=value
	elif opt in ("-e","--expression_fpkm_file"):
		expression_fpkm_file =value
	elif opt in ("-o","--out_dir"):
		out_dir =value  
	elif opt in ("-b","--binding_affinity_cutoff"):
		binding_affinity_cutoff =value
	elif opt in ("-a","--binding_affinity_foldchange_cutoff"):
		binding_affinity_foldchange_cutoff =value
	elif opt in ("-f","--fpkm_cutoff"):
		fpkm_cutoff =value 
	elif opt in ("-s","sample_id"):
		sample_id=value
	elif opt in ("-l","hla_str"):
		hla_str=value
#print coverage
if (input_netmhc_file =="" or input_fasta =="" or out_dir =="" or sample_id=="" or hla_str==""):
	print USAGE
	sys.exit(2)		
#######extract full animo acid change##
Full_header=[]
Full_gene=[]
Full_cdna=[]
Full_chrom_pos=[]
Full_transcript=[]
with open(input_fasta) as f:
	data=f.read()
for line in data.strip().split('\n'):
	if line.startswith('>WT'):
		full_aa_change=line.split('_')[2]
		full_gene=line.split('_')[1]
		full_cdna=line.split('_')[3]
		full_chr_p=line.split('_')[-2]
		full_transcript=line.split('_')[-1]
		Full_header.append(full_aa_change)
		Full_gene.append(full_gene)
		Full_cdna.append(full_cdna)
		Full_chrom_pos.append(full_chr_p)
		Full_transcript.append(full_transcript)
	else:
		continue
#print len(Full_header),len(Full_gene)
dup_full_header=[]
dup_full_gene=[]
dup_full_cdna=[]
dup_full_chrom_pos=[]
dup_full_transcript=[]
hla_num=len(hla_str.split(','))
i=0
while i<hla_num:
	for j in range(len(Full_header)):
		dup_full_header.append(Full_header[j])
		dup_full_gene.append(Full_gene[j])
		dup_full_cdna.append(Full_cdna[j])
		dup_full_chrom_pos.append(Full_chrom_pos[j])
		dup_full_transcript.append(Full_transcript[j])
	i=i+1

#print hla_num
#print Full_header
print len(dup_full_header),len(dup_full_gene)

######## extract candidate neoantigens####
#print expression_fpkm_file
with open(input_netmhc_file) as f:
    data = f.read()

nw_data = data.split('-----------------------------------------------------------------------------------\n')
#print len(nw_data)
#print len(nw_data)/8
WT_header = []
MT_header = []
WT_neo = []
MT_neo = []
for i in range(len(nw_data)):
    if i%8 == 3:
        wt_pro_name = nw_data[i].strip('\n').split('.')[0]
        WT_header.append(wt_pro_name)
    elif i%8 == 2:
        wt_neo_data = nw_data[i].strip().split('\n')
        WT_neo.append(wt_neo_data)
    elif i%8 == 7:
        mt_pro_name = nw_data[i].strip('\n').split('.')[0]
        MT_header.append(mt_pro_name)
    elif i%8 == 6:
        mt_neo_data = nw_data[i].strip().split('\n')
        MT_neo.append(mt_neo_data)
WB_SB_MT_record = []
WT_record = []
aa_record=[]
gene_record=[]
cdna_record=[]
chrom_pos_record=[]
transcript_name_record=[]
print len(MT_neo),len(WT_neo)


for i in range(len(MT_neo)):
	for j in range(len(MT_neo[i])):
		if MT_neo[i][j].endswith('WB') or MT_neo[i][j].endswith('SB'):
			#print i,j
			aa_record.append(dup_full_header[i])
			gene_record.append(dup_full_gene[i])
			#print dup_full_gene[i]
			cdna_record.append(dup_full_cdna[i])
			chrom_pos_record.append(dup_full_chrom_pos[i])
			transcript_name_record.append(dup_full_transcript[i])
			WB_SB_MT_record.append(MT_neo[i][j])
			WT_record.append(WT_neo[i][j])

candidate_neo = open(out_dir+'/'+sample_id+"_tmp_neo_candidate.txt",'w')
candidate_neo.write('\t'.join(['#Position','HLA_type','Gene','Transcript_name','Mutation','AA_change','MT_pep','WT_pep','MT_Binding_Aff','WT_Binding_Aff','MT_Binding_level_des','WT_Binding_level_des','fold_change']) + '\n')
for i in range(len(WB_SB_MT_record)):
    mt_record = [line for line in WB_SB_MT_record[i].split(' ') if line!='']
    HLA_tp = mt_record[1]
    gene = gene_record[i]
    transcript_n=transcript_name_record[i]
    ani_change = aa_record[i]
    cdna_change = cdna_record[i]
    chrom_pos_rec = chrom_pos_record[i]
    mt_pep = mt_record[2]
    mt_binding_rank=mt_record[13]
    mt_binding_level_des = mt_record[-1]
    wt_record = [i for i in WT_record[i].split(' ') if i!='']
    wt_pep = wt_record[2]
    wt_binding_rank=wt_record[13]
    if wt_record[-1]=='SB' or wt_record[-1]=='WB':
        wt_binding_level_des = wt_record[-1]
    else:
        wt_binding_level_des = 'NB'
    print wt_binding_rank
    print mt_binding_rank
    fold_change = float(wt_binding_rank)/float(mt_binding_rank)
    #DAI = float(wt_binding_aff) - float(mt_binding_aff)
    out_line = '\t'.join((chrom_pos_rec,HLA_tp,gene,transcript_n,cdna_change,ani_change,mt_pep,wt_pep,mt_binding_rank,wt_binding_rank,mt_binding_level_des,wt_binding_level_des,str(fold_change)))
    candidate_neo.write(out_line + '\n')
candidate_neo.close()
    

f=lambda x: x.split('.')[0]

######neoantigens filtering#####
##including Binding affinity ,localized peptide, multiple length peptide screen, differential AI> ##, neoantigens stability, gene FPKM>1 #####
data=pd.read_table(out_dir+'/'+sample_id+"_tmp_neo_candidate.txt",header=0,sep='\t')
if expression_fpkm_file=='no_exp':
	print "You did not provide expression file, the expression filter will not be done."
	final_filter_data=data[(data.MT_Binding_Aff<int(binding_affinity_cutoff))] 
elif os.path.exists(expression_fpkm_file):
	first_filter_data=data[(data.MT_Binding_Aff<int(binding_affinity_cutoff))]
	exp = pd.read_table(expression_fpkm_file,header=0,sep='\t')
	gene_exp = exp.loc[:,['target_id','tpm']]
	gene_exp["target_id"]=gene_exp.target_id.apply(f)
	neo_merge_exp = pd.merge(first_filter_data,gene_exp,left_on='Transcript_name',right_on='target_id',how='left')
	#print neo_merge_exp
	final_filter_data=neo_merge_exp[neo_merge_exp.tpm>float(fpkm_cutoff)]
else:
	print "could not find expresion file,check if the file exists!"
	sys.exit(2)
#os.remove(out_dir+'/'+sample_id+"_tmp_neo_candidate.txt")
#print final_filter_data
final_filter_data.to_csv(out_dir+'/'+sample_id+"_final_neo_candidate.txt",header=1,sep='\t',index=0)

    







            
