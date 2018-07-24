from pyfasta import Fasta
from collections import deque
import sys,getopt
opts,args=getopt.getopt(sys.argv[1:],"hi:c:o:s:",["input_eric_file","score_cutoff","out_dir","sample_id"])
input_eric_file=""
score_cutoff=0.5
out_dir=""
sample_id=""
USAGE='''
    This script convert insertion VCF derived VEP result to fasta format file for netMHC
    usage: python genefusion2fasta.py -i <input_vep_file> -c <score_cutoff> -o <outdir> -s <sample_id>
        required argument:
            -i | --input_eric_file : input file,result from ericscript
            -c | --score_cutoff : ericscript score cutoff for gene fusion detection
            -o | --out_dir : output directory
            -s | --sample_id : sample id
'''
for opt,value in opts:
    if opt =="h":
        print USAGE
        sys.exit(2)
    elif opt in ("-i","--input_eric_file"):
        input_eric_file=value
    elif opt in ("-c","--score_cutoff"):
        score_cutoff=value
    elif opt in ("-o","--out_dir"):
        out_dir=value
    elif opt in ("-s","--sample_id"):
        sample_id=value

#print coverage
if (input_eric_file =="" or out_dir =="" or sample_id==""):
    print USAGE
    sys.exit(2)

####
f_fasta=Fasta('/home/zhouchi/database/Annotation/Fasta/human.fasta')
global codon_dic
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

def reverse_nt(query_seq):
    d=deque()
    d.extendleft(query_seq)
    return ''.join(d)

header=[]
peptide=[]
with open(input_eric_file,'r') as f:
    data=f.read()

for long_line in data.strip().split('\n'):
    line=long_line.split('\t')
    head='>'+line[0]+'>>'+line[1]
    nt_seq_junction=line[18]
    ericscore=line[25]
    reverse_nt_seq_junction=reverse_nt(nt_seq_junction)
    codon_nt_forward_0=nt_seq_junction[21:81].upper()
    codon_nt_forward_1=nt_seq_junction[22:82].upper()
    codon_nt_forward_2=nt_seq_junction[23:83].upper()
    aa_forward_0=nt2aa(codon_nt_forward_0)
    aa_forward_1=nt2aa(codon_nt_forward_1)
    aa_forward_2=nt2aa(codon_nt_forward_2)
    #print aa_forward_0,aa_forward_1,aa_forward_2
    if ericscore >= score_cutoff:
        if len(aa_forward_0)>11:
            header.append(head)
            peptide.append(aa_forward_0)
        else:
            pass
        if len(aa_forward_1)>11:
            header.append(head)
            peptide.append(aa_forward_1)
        else:
            pass
        if len(aa_forward_2)>11:
            header.append(head)
            peptide.append(aa_forward_2)
        else:
            pass
    else: 
        continue

f_gf_o=open(out_dir+'/'+sample_id+'_gene_fusion.fasta','w')
for i in range(len(header)):
    f_gf_o.write('%s%s%s%s'%(header[i],'\n',peptide[i],'\n'))
f_gf_o.close()



