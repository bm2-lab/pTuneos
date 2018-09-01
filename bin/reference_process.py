f_o=open("human.fasta","w")

for line in open("Homo_sapiens_assembly38.fasta"):
	if line.startswith(">"):
		out=line.strip().split('  ')[0]
		f_o.write(out+'\n')
	else:
		f_o.write(line)

f_o.close()