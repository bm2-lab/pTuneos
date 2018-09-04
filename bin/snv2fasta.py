import pandas as pd
import numpy as np
import sys,getopt,os
import re
#####prepare fasta format input file for netMHC######
opts,args=getopt.getopt(sys.argv[1:],"hi:o:s:p:",["input_snv_vep_file","out_dir","sample_id","human_peptide"])
input_snv_vep_file =""
out_dir=""
sample_id=""
human_peptide=""
USAGE='''
	This script convert VEP result to fasta format file for netMHC
	usage: python animo_acid_prepare.py -i <input_snv_vep_file> -o <outdir> -s <sample_id> -p <human_peptide>
		required argument:
			-i | --input_snv_vep_file : input file,result from VEP
			-o | --out_dir : output directory
			-s | --sample_id : sample id
			-p | --human_peptide : reference protein sequence of human
'''
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-i","--input_snv_vep_file"):
		input_snv_vep_file=value
	elif opt in ("-o","--out_dir"):
		out_dir =value
	elif opt in ("-s","--sample_id"):
		sample_id =value 
	elif opt in ("-p","--human_peptide"):
		human_peptide =value 
	
#print coverage
if (input_snv_vep_file =="" or out_dir =="" or sample_id=="" or human_peptide==""):
	print USAGE
	sys.exit(2)		


transcript_aa={}
for line in open(human_peptide,'r'):
	if line.startswith(">"):
		transcript_name = line.strip().split(' ')[4][11:26]
		transcript_aa[transcript_name] = '' 
	else:
		transcript_aa[transcript_name] += line.replace('\n','')
protein_position=[]
cdna_position=[]
extra=[]
trans_name=[]
ref_animo_acid=[]
alt_animo_acid=[]
ref_nucleotide=[]
alt_nucleotide=[]
chrom_pos=[]
f_vep=open(input_snv_vep_file,'r')
for line in f_vep.readlines():
	if line.startswith('#'):
		pass
	else:
		record=line.strip().split('\t')
		chr_p=record[1]
		tran_n=record[4]
		pro_pos=record[9]
		ext=record[13]
		alt_aa=record[10].split('/')[1]
		ref_aa=record[10].split('/')[0]
		cdna_p=record[7]
		#print record[11].split('/')[0]
		ref_n=re.findall('[A-Z]',record[11].split('/')[0])[0]
		alt_n=re.findall('[A-Z]',record[11].split('/')[1])[0]
		alt_animo_acid.append(alt_aa)
		ref_animo_acid.append(ref_aa)
		trans_name.append(tran_n)
		protein_position.append(pro_pos)
		extra.append(ext)
		cdna_position.append(cdna_p)
		ref_nucleotide.append(ref_n)
		alt_nucleotide.append(alt_n)
		chrom_pos.append(chr_p)
f_vep.close()

gene_symbol=[]

for line in extra:
	element=line.split(';')
	for ele in element:
		sub_ele=ele.split('=')
		if sub_ele[0]=="SYMBOL":
			g_s= sub_ele[1]
	gene_symbol.append(g_s)

transcript_seq=[]
for name in trans_name:
	if name not in transcript_aa.keys():
		seq='NULL'	
	else:
		seq=transcript_aa[name]
	transcript_seq.append(seq)

mut_peptide=[]
wt_peptide=[]
mt_header=[]
wt_header=[]
for i in range(len(trans_name)):
	if transcript_seq[i]=="NULL":
		continue
	else:
		pro_change_pos=int(protein_position[i])
		wt_head='>WT_'+gene_symbol[i]+'_'+ ref_animo_acid[i]+str(pro_change_pos)+ alt_animo_acid[i]+'_'+ref_nucleotide[i]+str(cdna_position[i])+alt_nucleotide[i]+'_'+chrom_pos[i]+'_'+trans_name[i]
		mt_head='>MT_'+gene_symbol[i]+'_'+ ref_animo_acid[i]+str(pro_change_pos)+ alt_animo_acid[i]+'_'+ref_nucleotide[i]+str(cdna_position[i])+alt_nucleotide[i]+'_'+chrom_pos[i]+'_'+trans_name[i]
		ref_animo_acid_seq=transcript_seq[i]
		if pro_change_pos<=10:
			wt_pep=ref_animo_acid_seq[0:21]
			mt_pep=ref_animo_acid_seq[0:pro_change_pos-1]+alt_animo_acid[i]+ref_animo_acid_seq[pro_change_pos:21]
		elif pro_change_pos>10 and len(ref_animo_acid_seq)-pro_change_pos<=10:
			wt_pep=ref_animo_acid_seq[len(ref_animo_acid_seq)-21:len(ref_animo_acid_seq)]
			mt_pep=ref_animo_acid_seq[len(ref_animo_acid_seq)-21:pro_change_pos-1]+alt_animo_acid[i]+ref_animo_acid_seq[pro_change_pos:len(ref_animo_acid_seq)]
		else:
			wt_pep=ref_animo_acid_seq[pro_change_pos-11:pro_change_pos+10]
			mt_pep=ref_animo_acid_seq[pro_change_pos-11:pro_change_pos-1]+alt_animo_acid[i]+ref_animo_acid_seq[pro_change_pos:pro_change_pos+10]
		mt_header.append(mt_head)
		wt_header.append(wt_head)
		mut_peptide.append(mt_pep)
		wt_peptide.append(wt_pep)

mut_pep_len=[]
wt_pep_len=[]
for i in range(len(mut_peptide)):
	m_p_l=len(mut_peptide[i])
	w_p_l=len(wt_peptide[i])
	mut_pep_len.append(m_p_l)
	wt_pep_len.append(w_p_l)
#####drop duplicate###
snv_fasta_out=pd.DataFrame()
snv_fasta_out['mutation_header']=mt_header
snv_fasta_out['mutation_peptide']=mut_peptide
snv_fasta_out['wild_header']=wt_header
snv_fasta_out['wild_peptide']=wt_peptide
snv_fasta_out['mut_peptide_length']= mut_pep_len
snv_fasta_out['wt_peptide_length']= wt_pep_len
snv_fasta_dd=snv_fasta_out.drop_duplicates(subset=['mutation_header','mutation_peptide','wild_header','wild_peptide','mut_peptide_length'])
data_filter=snv_fasta_dd[(snv_fasta_dd["mut_peptide_length"]>=11) & (snv_fasta_dd["mut_peptide_length"]==snv_fasta_dd["wt_peptide_length"])]
data_snv_dd_reindex=data_filter.reset_index()
del data_snv_dd_reindex['index']
#######write######
f_w=open(out_dir+'/'+sample_id+"_snv.fasta",'w')
for i in range(len(data_snv_dd_reindex.mutation_header)):
	f_w.write('%s%s%s%s%s%s%s%s'%(data_snv_dd_reindex.wild_header[i],'\n',data_snv_dd_reindex.wild_peptide[i],'\n',data_snv_dd_reindex.mutation_header[i],'\n',data_snv_dd_reindex.mutation_peptide[i],'\n'))
f_w.close()











