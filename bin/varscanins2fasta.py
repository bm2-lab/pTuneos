from pyfasta import Fasta
import pandas as pd
import sys,getopt
opts,args=getopt.getopt(sys.argv[1:],"hi:o:s:r:p:",["input_ins_vep_file","out_dir","sample_id","reference","human_peptide_path"])
input_ins_vep_file =""
out_dir=""
sample_id=""
reference=""
human_peptide_path=""
USAGE='''
	This script convert insertion VCF derived VEP result to fasta format file for netMHC
	usage: python insertion2fasta.py -i <input_vep_file> -o <outdir> -s <sample_id> -r <reference> -p <human_peptide_path>
		required argument:
			-i | --input_ins_vep_file : input file,result from VEP
			-o | --out_dir : output directory
			-s | --sample_id : sample id
			-r | --reference : reference fasta file
			-p | --human_peptide_path : human_peptide_path
'''
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-i","--input_ins_vep_file"):
		input_ins_vep_file=value
	elif opt in ("-o","--out_dir"):
		out_dir=value
	elif opt in ("-s","--sample_id"):
		sample_id=value  
	elif opt in ("-r","--reference"):
		reference=value 
	elif opt in ("-p","--human_peptide_path"):
		human_peptide_path=value 	
#print coverage
if (input_ins_vep_file =="" or out_dir =="" or sample_id=="" or reference=="" or human_peptide_path==""):
	print USAGE
	sys.exit(2)	

####
f_fasta=Fasta(reference)
codon_dic={ 'TTT':'F','TTC':'F',
            'TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
            'TCT':'S','TCC':'S','TCA':'S','TCG':'S','AGT':'S','AGC':'S',
            'TAT':'Y','TAC':'Y',
            'TAA':'STOP','TAG':'STOP','TGA':'STOP',
            'TGT':'C','TGC':'C',
            'TGG':'W',
            'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
            'CAT':'H','CAC':'H',
            'CAA':'Q','CAG':'Q',
            'CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGA':'R','AGG':'R',
            'ATT':'I','ATC':'I','ATA':'I',
            'ATG':'M',
            'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
            'AAT':'N','AAC':'N',
            'AAA':'K','AAG':'K',
            'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
            'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
            'GAT':'D','GAC':'D',
            'GAA':'E','GAG':'E',
            'GGT':'G','GGC':'G','GGG':'G','GGA':'G'}

#############
transcript_aa={}
for line in open(human_peptide_path,'r'):
	if line.startswith(">"):
		transcript_name = line.strip().split(' ')[4][11:26]
		transcript_aa[transcript_name] = '' 
	else:
		transcript_aa[transcript_name] += line.replace('\n','')

def nt2aa(seq):
	aa_list=[]
	seq_len=len(seq)
	aa_num=seq_len/3
	for i in range(aa_num):
		nt_3=seq[i*3:i*3+3]
		nt_3_upper=nt_3.upper()
		if nt_3_upper in codon_dic.keys():
			aa_single=codon_dic[nt_3_upper]
			if aa_single!='STOP':
				aa_list.append(aa_single)
			else:
				break
		else:
			continue
	aa_seq=''.join(aa_list)
	return aa_seq
def nt_reverse(seq):
	seq_upper=seq.upper()
	nt_list=[]
	for i in range(len(seq_upper)):
		if seq_upper[i]=='A':
			nt_r='T'
		elif seq_upper[i]=='T':
			nt_r='A'
		elif seq_upper[i]=='C':
			nt_r='G'
		else:
			nt_r='C'
		nt_list.append(nt_r)
	nt_seq=''.join(nt_list)
	seq_reverse=nt_seq[::-1]
	return seq_reverse
location=[]
allele=[]
transcript_name=[]
consequence=[]
protein_position=[]
cds_position=[]
ref_animo_acid=[]
alt_animo_acid=[]
extra=[]
cdna_position=[]
animo_acid_change=[]
cdna_change=[]

for line in open(input_ins_vep_file):
	if line.startswith('#'):
		continue
	else:
		record=line.strip().split('\t')
		consequence_str=record[6].split(",")[0]
		if (consequence_str=="frameshift_variant" and record[0].split("_")[-1].split("/")[0]=="-") or (consequence_str=="inframe_insertion"):
			loc=record[1]
			alle=record[2]
			tran_n=record[4]
			cons=record[6].split(',')[0]
			pro_pos=record[9]
			cds_pos=record[8]
			if record[10]=='*':
				ref_aa='*'
				alt_aa='*'
			else:
				ref_aa=record[10].split('/')[0]
				alt_aa=record[10].split('/')[-1]
			ext=record[13]
			cdna_change_p=record[7]
			cdna_c=record[11]
			animo_acid_c=record[10]
			location.append(loc)
			allele.append(alle)
			transcript_name.append(tran_n)
			consequence.append(cons)
			protein_position.append(pro_pos)
			cds_position.append(cds_pos)
			ref_animo_acid.append(ref_aa)
			alt_animo_acid.append(alt_aa)
			extra.append(ext)
			animo_acid_change.append(animo_acid_c)
			cdna_position.append(cdna_change_p)
			cdna_change.append(cdna_c)
		else:
			continue
gene_symbol=[]
strand=[]
for line in extra:
	element=line.split(';')
	for ele in element:
		sub_ele=ele.split('=')
		if sub_ele[0]=="SYMBOL":
			g_s=sub_ele[1]
		elif sub_ele[0]=="STRAND":
			st=sub_ele[1]
	gene_symbol.append(g_s)
	strand.append(st)
transcript_seq=[]
for name in transcript_name:
	if name not in transcript_aa.keys():
		seq='NULL'
	else:
		seq=transcript_aa[name]
	transcript_seq.append(seq)
mut_piptide=[]
wt_piptide=[]
mt_header=[]
wt_header=[]
for i in range(len(location)):
	if transcript_seq[i]=="NULL":
		continue
	elif consequence[i]=='inframe_insertion':
		chr_name=location[i].split(':')[0]
		ins_start=int(location[i].split(':')[1].split('-')[0])
		cds_start=int(cds_position[i].split('-')[0])
		protein_change_pos=int(protein_position[i].split('-')[0])
		ref_aa=ref_animo_acid[i]
		alt_aa=alt_animo_acid[i]
		trans_strand=int(strand[i])
		seq=transcript_aa[transcript_name[i]]
		wt_head='>WT_'+gene_symbol[i]+'_'+protein_position[i]+'.'+animo_acid_change[i]+'_'+cdna_position[i]+'.'+cdna_change[i]+'_'+location[i]+'_'+transcript_name[i]
		mt_head='>MT_'+gene_symbol[i]+'_'+protein_position[i]+'.'+animo_acid_change[i]+'_'+cdna_position[i]+'.'+cdna_change[i]+'_'+location[i]+'_'+transcript_name[i]
		alt_aa_len=len(alt_aa)
		if int(protein_change_pos)<=10:
			wt_pt=seq[0:protein_change_pos+10]
			if alt_aa_len>=11:
				mt_pt=seq[0:protein_change_pos-1]+alt_aa[0:11]
			else:
				mt_pt=seq[0:protein_change_pos-1]+alt_aa+seq[protein_change_pos:protein_change_pos+11-len(alt_aa)]
		elif protein_change_pos>10 and len(seq)-protein_change_pos<=10:
			wt_pt=seq[protein_change_pos-11:len(seq)]
			if alt_aa_len>=11:
				mt_pt=seq[protein_change_pos-11:protein_change_pos-1]+alt_aa[0:len(seq)-protein_change_pos]
			elif alt_aa_len>=(len(seq)-protein_change_pos):
				mt_pt=seq[protein_change_pos-11:protein_change_pos-1]+alt_aa[0:len(seq)-protein_change_pos]
			else:
				mt_pt=seq[protein_change_pos-11:protein_change_pos-1]+alt_aa+seq[len(seq)-protein_change_pos-alt_aa_len:len(seq)]
		else:
			wt_pt=seq[protein_change_pos-11:protein_change_pos+10]
			if ref_aa=='-':
				if alt_aa_len>=11:
					mt_pt=seq[protein_change_pos-11:protein_change_pos-1]+alt_aa[0:11]
				else:
					mt_pt=seq[protein_change_pos-11:protein_change_pos-1]+alt_aa+seq[protein_change_pos-1:protein_change_pos+10-len(alt_aa)]
			else:
				if alt_aa_len>=11:
					mt_pt=seq[protein_change_pos-11:protein_change_pos-1]+alt_aa[0:11]
				else:
					mt_pt=seq[protein_change_pos-11:protein_change_pos-1]+alt_aa+seq[protein_change_pos:protein_change_pos+11-len(alt_aa)]
		mut_piptide.append(mt_pt)
		wt_piptide.append(wt_pt)
		mt_header.append(mt_head)
		wt_header.append(wt_head)
	elif consequence[i]=='frameshift_variant':
		chr_name=location[i].split(':')[0]
		ins_start=int(location[i].split(':')[1].split('-')[0])
		cds_start=int(cds_position[i].split('-')[0])
		protein_change_pos=int(protein_position[i].split('-')[0])
		allele_ins=allele[i]
		trans_strand=int(strand[i])
		frame_left_num=int(cds_start)%3
		seq=transcript_aa[transcript_name[i]]
		if protein_change_pos<11:
			aa_ten_before=seq[0:protein_change_pos-1]
		else:
			aa_ten_before=seq[protein_change_pos-11:protein_change_pos-1]
		wt_head='>WT_'+gene_symbol[i]+'_'+protein_position[i]+'.'+animo_acid_change[i]+'_'+cdna_position[i]+'.'+cdna_change[i]+'_'+location[i]+'_'+transcript_name[i]
		mt_head='>MT_'+gene_symbol[i]+'_'+protein_position[i]+'.'+animo_acid_change[i]+'_'+cdna_position[i]+'.'+cdna_change[i]+'_'+location[i]+'_'+transcript_name[i]
		#print wt_head
		if trans_strand==1:
			if len(allele_ins)>=33-frame_left_num:
				nt_all=f_fasta[chr_name][ins_start-frame_left_num:ins_start]+allele_ins[0:33-frame_left_num]
				#print nt_all
				mt_pt_all=aa_ten_before+nt2aa(nt_all)
				#print mt_pt
				wt_pt=seq[protein_change_pos-11:protein_change_pos-1+len(nt2aa(nt_all))]
				if len(mt_pt_all)>len(wt_pt):
					mt_pt=mt_pt_all[0:len(wt_pt)]
				else:
					mt_pt=mt_pt_all
				#wt_nt=f_fasta[chr_name][ins_start-31:ins_start+32]
				#wt_pt=nt2aa(wt_nt)
				#print wt_nt
				#print wt_pt
				#print mt_pt
			else:
				nt_all=f_fasta[chr_name][ins_start-frame_left_num:ins_start]+allele_ins+f_fasta[chr_name][ins_start:ins_start+33-frame_left_num-len(allele_ins)]
				#print len(nt_all)
				#print aa_ten_before,nt2aa(nt_all)
				mt_pt_all=aa_ten_before+nt2aa(nt_all)
				wt_pt=seq[protein_change_pos-11:protein_change_pos-1+len(nt2aa(nt_all))]
				if len(mt_pt_all)>len(wt_pt):
					mt_pt=mt_pt_all[0:len(wt_pt)]
				else:
					mt_pt=mt_pt_all
		else:
			reverse_ins=nt_reverse(allele_ins)
			if len(allele_ins)>=33-frame_left_num:
				nt_all=nt_reverse(f_fasta[chr_name][ins_start:ins_start+frame_left_num])+reverse_ins[0:33-frame_left_num]
				#print nt_all
				mt_pt_all=aa_ten_before+nt2aa(nt_all)
				wt_pt=seq[protein_change_pos-11:protein_change_pos-1+len(nt2aa(nt_all))]
				if len(mt_pt_all)>len(wt_pt):
					mt_pt=mt_pt_all[0:len(wt_pt)]
				else:
					mt_pt=mt_pt_all

			else:
				nt_all=nt_reverse(f_fasta[chr_name][ins_start:ins_start+frame_left_num])+reverse_ins+nt_reverse(f_fasta[chr_name][ins_start+frame_left_num+len(reverse_ins)-34:ins_start-1])
				#print len(nt_all),nt_all,nt2aa(nt_all)
				mt_pt_all=aa_ten_before+nt2aa(nt_all)
				wt_pt=seq[protein_change_pos-11:protein_change_pos-1+len(nt2aa(nt_all))]
				if len(mt_pt_all)>len(wt_pt):
					mt_pt=mt_pt_all[0:len(wt_pt)]
				else:
					mt_pt=mt_pt_all
		if mt_pt!='' and wt_pt!='':
			mut_piptide.append(mt_pt)
			wt_piptide.append(wt_pt)
			mt_header.append(mt_head)
			wt_header.append(wt_head)
		else:
			continue
	else:
		continue

###drop duplicate###

mut_pep_len=[]
wt_pep_len=[]
for i in range(len(mut_piptide)):
	m_p_l=len(mut_piptide[i])
	mut_pep_len.append(m_p_l)
	w_p_l=len(wt_piptide[i])
	wt_pep_len.append(w_p_l)

ins_fasta_out=pd.DataFrame()
ins_fasta_out['mutation_header']=mt_header
ins_fasta_out['mutation_peptide']=mut_piptide
ins_fasta_out['wild_header']=wt_header
ins_fasta_out['wild_peptide']=wt_piptide
ins_fasta_out['mutation_peptide_length']=mut_pep_len
ins_fasta_out['wild_peptide_length']=wt_pep_len
ins_fasta_dd=ins_fasta_out.drop_duplicates(subset=['mutation_header','mutation_peptide','wild_header','wild_peptide'])
data_filter=ins_fasta_dd[(ins_fasta_dd["mutation_peptide_length"]>=11) & (ins_fasta_dd["mutation_peptide_length"]==ins_fasta_dd["wild_peptide_length"])]
#data_filter=ins_fasta_dd
#data_filter.dropna(axis=0)
data_ins_dd_reindex=data_filter.reset_index()
del data_ins_dd_reindex['index']
f_w=open(out_dir+'/'+sample_id+'_ins.fasta','w')
for i in range(len(data_ins_dd_reindex.mutation_header)):
	f_w.write('%s%s%s%s%s%s%s%s'%(data_ins_dd_reindex.wild_header[i],'\n',data_ins_dd_reindex.wild_peptide[i],'\n',data_ins_dd_reindex.mutation_header[i],'\n',data_ins_dd_reindex.mutation_peptide[i],'\n'))
f_w.close()







