import os

f_o=open("human.fasta","w")
if os.path.exists("Homo_sapiens_assembly38.fasta"):
	for line in open("Homo_sapiens_assembly38.fasta"):
		if line.startswith(">"):
			out=line.strip().split('  ')[0]
			f_o.write(out+'\n')
		else:
			f_o.write(line)
else:
	for line in open("ucsc.hg19.fasta"):
		if line.startswith(">"):
			out=line.strip().split('  ')[0]
			f_o.write(out+'\n')
		else:
			f_o.write(line)		

f_o.close()