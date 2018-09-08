from pyfasta import Fasta
import pandas as pd
import sys,getopt
opts,args=getopt.getopt(sys.argv[1:],"hi:o:s:",["input_del_vep_file","out_dir","sample_id"])
input_del_vep_file =""
out_dir=""
sample_id=""
USAGE='''
	This script convert deletion VCF derived VEP result to fasta format file for netMHC
	usage: python deletion2fasta.py -i <input_vep_file> -o <outdir> -s <sample_id>
		required argument:
			-i | --input_del_vep_file : input file,result from VEP
			-o | --out_dir : output directory
			-s | --sample_id : sample id
'''
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-i","--input_del_vep_file"):
		input_del_vep_file=value
	elif opt in ("-o","--out_dir"):
		out_dir =value
	elif opt in ("-s","--sample_id"):
		sample_id =value  
	
#print coverage
if (input_del_vep_file =="" or out_dir =="" or sample_id==""):
	print USAGE
	sys.exit(2)	
####
f_fasta=Fasta('/home/zhouchi/database/Annotation/Fasta/human.fasta')
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
for line in open("/home/zhouchi/database/Annotation/protein/Homo_sapiens.GRCh38.pep.all.fa",'r'):
	if line.startswith(">"):
		transcript_name = line.strip().split(' ')[4][11:26]
		transcript_aa[transcript_name] = '' 
	else:
		transcript_aa[transcript_name] += line.replace('\n','')

location=[]
allele=[]
transcript_name=[]
consequence=[]
protein_position=[]
cds_position=[]
extra=[]
cdna_position=[]
animo_acid_change=[]
cdna_change=[]
f_vep=open(input_del_vep_file,'r')
for line in f_vep.readlines():
	if line.startswith('#'):
		pass
	else:
		record=line.strip().split('\t')
		loc=record[1]
		alle=record[2]
		tran_n=record[4]
		cons=record[6].split(',')[0]
		pro_pos=record[9]
		cds_pos=record[8]
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
		extra.append(ext)
		animo_acid_change.append(animo_acid_c)
		cdna_position.append(cdna_change_p)
		cdna_change.append(cdna_c)
f_vep.close()
gene_symbol=[]
strand=[]
for line in extra:
	element=line.split(';')
	for ele in element:
		sub_ele=ele.split('=')
		if sub_ele[0]=="SYMBOL":
			g_s=sub_ele[1]
		if sub_ele[0]=="STRAND":
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
def nt_reverse(seq):
	reverse_list=[]
	seq_f_list=list(seq)
	for ele in seq_f_list:
		if ele=='A':
			reverse_ele='T'
		elif ele=='T':
			reverse_ele='A'
		elif ele=='C':
			reverse_ele='G'
		else:
			reverse_ele='C'
		reverse_list.append(reverse_ele)
	seq_str=''.join(reverse_list)
	reverse_seq=seq_str[::-1]
	return reverse_seq
def nt2aa(seq):
	aa_list=[]
	seq_len=len(seq)
	aa_num=seq_len/3
	for i in range(aa_num):
		nt_3=seq[i*3:i*3+3]
		nt_3_upper=nt_3.upper()
		aa_single=codon_dic[nt_3_upper]
		if aa_single!='STOP':
			aa_list.append(aa_single)
		else:
			break
	aa_seq=''.join(aa_list)
	return aa_seq
for i in range(len(location)):
	if transcript_seq[i]=="NULL":
		continue
	else:
		if consequence[i]=="inframe_deletion":
			chr_name=location[i].split(':')[0]
			del_start=location[i].split(':')[1].split('-')[0]
			del_end=location[i].split(':')[1].split('-')[-1]
			cds_start=cds_position[i].split('-')[0]
			cds_end=cds_position[i].split('-')[-1]
			protein_change_pos= protein_position[i].split('-')[0]
			strand_n=int(strand[i])
		#	protein_change_pos_end= protein_position[i].split('-')[-1]
			frame_left_num=(int(cds_start)-1)%3
			seq=transcript_aa[transcript_name[i]]
			wt_head='>WT_'+gene_symbol[i]+'_'+protein_position[i]+'.'+animo_acid_change[i]+'_'+cdna_position[i]+'.'+cdna_change[i]+'_'+location[i]+'_'+transcript_name[i]
			mt_head='>MT_'+gene_symbol[i]+'_'+protein_position[i]+'.'+animo_acid_change[i]+'_'+cdna_position[i]+'.'+cdna_change[i]+'_'+location[i]+'_'+transcript_name[i]
			if frame_left_num==0:
				del_pro_num=(int(del_end)-int(del_start)+1)/3
				if int(protein_change_pos)<=10:
					mt_pt=seq[0:(int(protein_change_pos)-1)]+seq[int(protein_change_pos)+del_pro_num:int(protein_change_pos)+del_pro_num+(22-int(protein_change_pos))]
					wt_pt=seq[0:21]
				elif (int(protein_change_pos)>10 and len(seq)-int(protein_change_pos)<=10):
					mt_pt=seq[len(seq)-21-del_pro_num:(int(protein_change_pos)-1)]+seq[int(protein_change_pos)+del_pro_num:len(seq)+1]
					wt_pt=seq[len(seq)-21:len(seq)+1]
				else:
					mt_pt=seq[int(protein_change_pos)-11:int(protein_change_pos)-1]+seq[int(protein_change_pos)+del_pro_num-1:int(protein_change_pos)+del_pro_num+10]
					wt_pt=seq[int(protein_change_pos)-11:int(protein_change_pos)+10]
				mut_piptide.append(mt_pt)
				wt_piptide.append(wt_pt)
				mt_header.append(mt_head)
				wt_header.append(wt_head)
			elif frame_left_num==1:
				del_pro_num=(int(del_end)-int(del_start)+1)/3
				if strand_n==1:
					change_nt_l=f_fasta[chr_name][int(del_start)-2]
					change_nt_r=f_fasta[chr_name][int(del_end):int(del_end)+2]
					change_nt=change_nt_l+change_nt_r
					change_nt_upper=change_nt.upper()
					change_aa=codon_dic[change_nt_upper]
				else:
					change_nt_l=f_fasta[chr_name][int(del_start)-3:int(del_start)-1]
					change_nt_r=f_fasta[chr_name][int(del_end):int(del_end)+1]
					change_nt=change_nt_l+change_nt_r
					change_nt_upper=change_nt.upper()
					change_nt_upper_reverse=nt_reverse(change_nt_upper)
					change_aa=codon_dic[change_nt_upper_reverse]					
				if change_aa!='STOP':
					if int(protein_change_pos)<=10:
						mt_pt=seq[0:int(protein_change_pos)-1]+change_aa+seq[int(protein_change_pos)+del_pro_num:int(protein_change_pos)+del_pro_num+21-int(protein_change_pos)]
						wt_pt=seq[0:21]
					elif (int(protein_change_pos)>10 and len(seq)-int(protein_change_pos)<=10):
						mt_pt=seq[len(seq)-21-del_pro_num:(int(protein_change_pos)-1)]+seq[int(protein_change_pos)+del_pro_num:len(seq)]
						wt_pt=wt_pt=seq[len(seq)-21:len(seq)+1]
					else:
						mt_pt=seq[int(protein_change_pos)-11:int(protein_change_pos)-1]+change_aa+seq[int(protein_change_pos)+del_pro_num:int(protein_change_pos)+del_pro_num+10]
						wt_pt=seq[int(protein_change_pos)-11:int(protein_change_pos)+10]
					mut_piptide.append(mt_pt)
					wt_piptide.append(wt_pt)
					mt_header.append(mt_head)
					wt_header.append(wt_head)
				else:
					print '1the deletion result in a STOP codon,so no peptide generate!'
			elif frame_left_num==2:
				del_pro_num=(int(del_end)-int(del_start)+1)/3
				if strand_n==1:
					change_nt_l=f_fasta[chr_name][int(del_start)-3:int(del_start)-1]
					change_nt_r=f_fasta[chr_name][int(del_end)]
					change_nt=change_nt_l+change_nt_r
					change_nt_upper=change_nt.upper()
					change_aa=codon_dic[change_nt_upper]
				else:
					change_nt_l=f_fasta[chr_name][int(del_start)-2]
					change_nt_r=f_fasta[chr_name][int(del_end):int(del_end)+2]
					change_nt=change_nt_l+change_nt_r
					change_nt_upper=change_nt.upper()
					change_nt_upper_reverse=nt_reverse(change_nt_upper)
					change_aa=codon_dic[change_nt_upper_reverse]
				if change_aa!='STOP':
					if int(protein_change_pos)<=10:
						mt_pt=seq[0:int(protein_change_pos)-1]+change_aa+seq[int(protein_change_pos)+del_pro_num:int(protein_change_pos)+del_pro_num+21-int(protein_change_pos)]
						wt_pt=seq[0:21]
					elif (int(protein_change_pos)>10 and len(seq)-int(protein_change_pos)<=10):
						mt_pt=seq[len(seq)-21-del_pro_num:(int(protein_change_pos)-1)]+seq[int(protein_change_pos)+del_pro_num:len(seq)]
						wt_pt=seq[len(seq)-21:len(seq)+1]
					else:
						mt_pt=seq[int(protein_change_pos)-11:int(protein_change_pos)-1]+change_aa+seq[int(protein_change_pos)+del_pro_num:int(protein_change_pos)+del_pro_num+10]
						wt_pt=seq[int(protein_change_pos)-11:int(protein_change_pos)+10]
					mut_piptide.append(mt_pt)
					wt_piptide.append(wt_pt)
					mt_header.append(mt_head)
					wt_header.append(wt_head)
				else:
					print '2the deletion result in a STOP codon,so no peptide generate!'
				#mut_piptide.append(mt_pt)
				#wt_piptide.append(wt_pt)
		if consequence[i]=="frameshift_variant":
			chr_name=location[i].split(':')[0]
			del_loc_start=int(location[i].split(':')[1].split('-')[0])
			del_loc_end=int(location[i].split(':')[1].split('-')[-1])
			#print del_loc_start
			strand_n=int(strand[i])
			cds_loc=cds_position[i].split('-')[0]
			#print protein_position[i],consequence[i]
			protein_change_pos_start=int(protein_position[i].split('-')[0])
			frame_left_num=(int(cds_loc)-1)%3
			seq=transcript_aa[transcript_name[i]]
			wt_head='>WT_'+gene_symbol[i]+'_'+protein_position[i]+'.'+animo_acid_change[i]+'_'+cdna_position[i]+'.'+cdna_change[i]+'_'+location[i]+'_'+transcript_name[i]
			mt_head='>MT_'+gene_symbol[i]+'_'+protein_position[i]+'.'+animo_acid_change[i]+'_'+cdna_position[i]+'.'+cdna_change[i]+'_'+location[i]+'_'+transcript_name[i]
			if strand_n==1:
				nt_left=f_fasta[chr_name][del_loc_start-1-frame_left_num:del_loc_start-1]
				if int(protein_change_pos_start)<=10:
					nt_right=f_fasta[chr_name][del_loc_end:del_loc_end+34-frame_left_num]
					nt_change=nt_left+nt_right
					aa_change=nt2aa(nt_change)
					mt_pt=seq[0:protein_change_pos_start-1]+aa_change
					mt_aa_len=len(mt_pt)
					wt_pt=seq[0:mt_aa_len]
				elif (int(protein_change_pos_start)>10 and len(seq)-int(protein_change_pos_start)<=10):
					right_aa_num=len(seq)-protein_change_pos_start
					nt_right=f_fasta[chr_name][del_loc_end+1:del_loc_end+right_aa_num*3]
					nt_change=nt_left+nt_right
					aa_change=nt2aa(nt_change)
					mt_pt=seq[int(protein_change_pos_start)-11:protein_change_pos_start-1]+aa_change
					wt_pt=seq[int(protein_change_pos_start)-11:protein_change_pos_start-1+len(aa_change)]
				else:
					nt_right=f_fasta[chr_name][del_loc_end:del_loc_end+33-frame_left_num]
					nt_change=nt_left+nt_right
					#print nt_change
					#print len(nt_change)
					aa_change=nt2aa(nt_change)
					#print len(aa_change),aa_change
					mt_pt=seq[int(protein_change_pos_start)-11:protein_change_pos_start-1]+aa_change
					mt_aa_len=len(aa_change)
					wt_pt=seq[int(protein_change_pos_start)-11:int(protein_change_pos_start)-1+mt_aa_len]
				mut_piptide.append(mt_pt)
				wt_piptide.append(wt_pt)
				mt_header.append(mt_head)
				wt_header.append(wt_head)
			else:
				change_nt_right_all_reverse=nt_reverse(f_fasta[chr_name][del_loc_start-(3-frame_left_num)-31:del_loc_start+frame_left_num])
				#print change_nt_right_all_reverse
				if frame_left_num==0:
					change_nt=change_nt_right_all_reverse[1:]
				elif frame_left_num==1:
					change_nt=change_nt_right_all_reverse[0]+change_nt_right_all_reverse[2:]
				else:
					change_nt=change_nt_right_all_reverse[0:2]+change_nt_right_all_reverse[3:]
				#print change_nt
				if int(protein_change_pos_start)<=10:
					aa_change=nt2aa(change_nt)
					mt_pt=seq[0:protein_change_pos_start-1]+aa_change
					mt_aa_len=len(mt_pt)
					wt_pt=seq[0:mt_aa_len]
				elif (int(protein_change_pos_start)>10 and len(seq)-int(protein_change_pos_start)<=10):
					aa_change=nt2aa(change_nt)
					mt_pt=seq[int(protein_change_pos_start)-11:protein_change_pos_start]+aa_change
					wt_pt=seq[int(protein_change_pos_start)-11:protein_change_pos_start+len(aa_change)]
				else:
					aa_change=nt2aa(change_nt)
					#print aa_change
					mt_pt=seq[int(protein_change_pos_start)-11:protein_change_pos_start-1]+aa_change
					mt_aa_len=len(aa_change)
					wt_pt=seq[int(protein_change_pos_start)-11:int(protein_change_pos_start)-1+mt_aa_len]
				mut_piptide.append(mt_pt)
				wt_piptide.append(wt_pt)
				mt_header.append(mt_head)
				wt_header.append(wt_head)
		else:
			pass


mut_pep_len=[]
wt_pep_len=[]
for i in range(len(mut_piptide)):
	m_p_l=len(mut_piptide[i])
	w_p_l=len(wt_piptide[i])
	mut_pep_len.append(m_p_l)
	wt_pep_len.append(w_p_l)
###drop duplicate###
del_fasta_out=pd.DataFrame()
del_fasta_out['mutation_header']=mt_header
del_fasta_out['mutation_peptide']=mut_piptide
del_fasta_out['wild_header']=wt_header
del_fasta_out['wild_peptide']=wt_piptide
del_fasta_out['mut_peptide_length']= mut_pep_len
del_fasta_out['wt_peptide_length']= wt_pep_len
del_fasta_dd=del_fasta_out.drop_duplicates(subset=['mutation_header','mutation_peptide','wild_header','wild_peptide','mut_peptide_length'])
data_filter=del_fasta_dd[(del_fasta_dd["mut_peptide_length"]>=11) & (del_fasta_dd["mut_peptide_length"]==del_fasta_dd["wt_peptide_length"])]
data_del_dd_reindex=data_filter.reset_index()
del data_del_dd_reindex['index']
#######write######
f_w=open(out_dir+'/'+sample_id+"_del.fasta",'w')
for i in range(len(data_del_dd_reindex.mutation_header)):
	f_w.write('%s%s%s%s%s%s%s%s'%(data_del_dd_reindex.wild_header[i],'\n',data_del_dd_reindex.wild_peptide[i],'\n',data_del_dd_reindex.mutation_header[i],'\n',data_del_dd_reindex.mutation_peptide[i],'\n'))
f_w.close()













