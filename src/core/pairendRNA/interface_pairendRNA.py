import os
from pairendRNAprocessor import *
import shutil
import yaml
def PERNA(opts):
	config_file=opts.Config_file
	f=open(config_file)
	config_list=yaml.load(f)
	#######read and parse parameter
	print "read and parse parameter."
	output_fold=config_list["output_fold"]
	itunes_bin_path=config_list["itunes_bin_path"]
	os.system("export iTuNES_BIN_PATH=%s"%itunes_bin_path)
	tumor_fastq_path_first=config_list["tumor_fastq_path_first"]
	tumor_fastq_path_second=config_list["tumor_fastq_path_second"]
	opitype_fold=config_list["opitype_fold"]
	opitype_out_fold=output_fold + '/' + 'hlatyping'
	opitype_ext='${iTuNES_BIN_PATH}/optitype_ext.py'
	prefix=config_list["sample_name"]
	CPU=config_list["thread_number"]
	vcftools_path=config_list["vcftools_path"]
	REFERENCE=config_list["reference_path"]
	STAR_INDEX=config_list["STAR_index_path"]
	GENOME=config_list["genome"]
	alignment_out_fold=output_fold + '/' + 'alignments'
	STAR_path=config_list["STAR_path"]
	samtools_path=config_list["samtools_path"]
	java_picard_path=config_list["java_picard_path"]
	GATK_path=config_list["GATK_path"]
	dbsnp138_path=config_list["dbsnp138"]
	hapmap_path=config_list["hapmap"]
	omni_path=config_list["omni"]
	OneKG_path=config_list["OneKG"]
	mills_path=config_list["mills"]
	cosmic_path=config_list["cosmic"]
	somatic_mutation_fold=output_fold + '/' + 'somatic_mutation'
	vep_cache=config_list["vep_cache_path"]
	varscan_path=config_list["varscan_path"]
	vep_path=config_list["vep_path"]
	netmhc_out_fold=output_fold + '/' + 'netmhc'
	varscan_indel_fold=output_fold + '/' + 'indel'
	strelka_out_fold=output_fold + '/' + 'strelka'
	indel_fasta_file=netmhc_out_fold+'/'+prefix+'_indel.fasta'
	hla_str=config_list["hla_str"]
	indel_netmhc_out_file=netmhc_out_fold+'/'+prefix+'_indel_netmhc.txt'
	split_num=500
	binding_fc_aff_cutoff=int(config_list["binding_fc_aff_cutoff"])
	binding_aff_cutoff=int(config_list["binding_aff_cutoff"])
	fpkm_cutoff=int(config_list["fpkm_cutoff"])
	netctl_out_fold=output_fold + '/' + 'netctl'
	netMHCpan_path=config_list["netMHCpan_path"]
	snv_fasta_file=netmhc_out_fold+'/'+prefix+'_snv.fasta'
	snv_netmhc_out_file=netmhc_out_fold+'/'+prefix+'_snv_netmhc.txt'
	strelka_path=config_list["strelka_path"]
	kallisto_path=config_list["kallisto_path"]
	kallisto_cdna_path=config_list["kallisto_cdna_path"]
	final_neo_file=netctl_out_fold + '/' + prefix + '_snv_netctl_concact.txt'
	gmm_classification_file=netctl_out_fold + '/' + prefix + '_gmm.png'
	immunogenicity_score_ranking=netctl_out_fold + '/' + prefix + '_score.txt'
	immunogenicity_gmm_score_ranking=netctl_out_fold + '/' + prefix + '_gmm_score.txt'
	candidate_neoantigens_fold=output_fold + '/' + 'candidate_neoantigens'
	GTF_path=config_list["GTF_path"]
	kallisto_out_fold=output_fold + '/' + 'expression'
	logfile_out_fold=output_fold + '/' + 'logfile'
	clean_fastq_fold=output_fold + '/' + 'clean_fastq'
	trimmomatic_path=config_list["trimmomatic_path"]
	adapter_path=config_list["adapter_path"]
	tumor_fastq_prefix=clean_fastq_fold + '/' + prefix + "_tumor_clean.fq.gz"
	normal_fastq_prefix=clean_fastq_fold + '/' + prefix + "_normal_clean.fq.gz"
	tumor_fastq_clean_first=clean_fastq_fold + '/' + prefix + "_tumor_clean_1P.fq.gz"
	tumor_fastq_clean_second=clean_fastq_fold + '/' + prefix + "_tumor_clean_2P.fq.gz"
	pos_1000G_file_path=config_list["pos_1000G_file"]
	#####check input file,tool path and reference file#####
	if os.path.exists(tumor_fastq_path_first) and os.path.exists(tumor_fastq_path_second):
		print "check all fastq file done."
	else:
		print "please check your input fastq file!"
		os._exit(1)
	if os.path.exists(opitype_fold):
		print "check opitype path done."
	else:
		print "please check your opitype path!"
		os._exit(1)
	if os.path.exists(STAR_path):
		print "check STAR path done."
	else:
		print "please check your STAR path!"
		os._exit(1)	
	if os.path.exists(java_picard_path):
		print "check picard path done."
	else:
		print "please check your picard path!"
		os._exit(1)	
	if os.path.exists(GATK_path):
		print "check GATK path done."
	else:
		print "please check your GATK path!"
		os._exit(1)	
	if os.path.exists(dbsnp138_path):
		print "check dbsnp138 file path done."
	else:
		print "please check your dbsnp138 file path!"
		os._exit(1)	
	if os.path.exists(OneKG_path):
		print "check OneKG file path done."
	else:
		print "please check your OneKG file path!"
		os._exit(1)	
	if os.path.exists(mills_path):
		print "check mills file path done."
	else:
		print "please check your mills file path!"
		os._exit(1)	
	if os.path.exists(vep_path):
		print "check vep path done."
	else:
		print "please check your vep path!"
		os._exit(1)	
	if os.path.exists(vep_cache):
		print "check vep cache path done."
	else:
		print "please check your vep cache path!"
		os._exit(1)	
	exp_file=output_fold + '/' + "expression/abundance.tsv"
	if os.path.exists(REFERENCE):
		print "check REFERENCE file path done."
	else:
		print "please check your REFERENCE file path!"
		os._exit(1)	
	if os.path.exists(STAR_INDEX+'/'+"transcriptInfo.tab") and os.path.exists(STAR_INDEX+'/'+"sjdbList.out.tab") and os.path.exists(STAR_INDEX+'/'+"sjdbList.fromGTF.out.tab") and os.path.exists(STAR_INDEX+'/'+"sjdbInfo.txt") and os.path.exists(STAR_INDEX+'/'+"SAindex") and os.path.exists(STAR_INDEX+'/'+"SA") and os.path.exists(STAR_INDEX+'/'+"genomeParameters.txt") and os.path.exists(STAR_INDEX+'/'+"Genome") and os.path.exists(STAR_INDEX+'/'+"geneInfo.tab"): 
		print "STAR index path:%s"%(STAR_INDEX)
	elif os.path.exists(REFERENCE):
		print "no bwa index file, make index."
		if not os.path.exists(STAR_INDEX):
			os.mkdir(STAR_INDEX)
		cmd_make_index=STAR_path+" --runMode genomeGenerate --runThreadN 8 --genomeDir "+STAR_INDEX + " --genomeFastaFiles "+ REFERENCE + " --sjdbGTFfile " + GTF_path + " --sjdbOverhang 99"
		print cmd_make_index
		os.system(cmd_make_index)
	else:
		print "ERROR: no index file and no genome.fa to build it"
		os._exit(1)
	#####check output directory###
	print "check output directory"
	if not os.path.exists(output_fold):
		os.mkdir(output_fold)
	if not os.path.exists(somatic_mutation_fold):
		os.mkdir(somatic_mutation_fold)
	if not os.path.exists(netmhc_out_fold):
		os.mkdir(netmhc_out_fold)
	if not os.path.exists(netctl_out_fold):
		os.mkdir(netctl_out_fold)
	if not os.path.exists(alignment_out_fold):
		os.mkdir(alignment_out_fold)
	if not os.path.exists(kallisto_out_fold):
		os.mkdir(kallisto_out_fold)
	if not os.path.exists(clean_fastq_fold):
		os.mkdir(clean_fastq_fold)
	if not os.path.exists(candidate_neoantigens_fold):
		os.mkdir(candidate_neoantigens_fold)
	if not os.path.exists(candidate_neoantigens_fold):
		os.mkdir(candidate_neoantigens_fold)
	if not os.path.exists(logfile_out_fold):
		os.mkdir(logfile_out_fold)
	if not os.path.exists(opitype_out_fold):
		os.mkdir(opitype_out_fold)
	print "start fastq quality control"
	processes_0=[]
	q1=multiprocessing.Process(target=read_trimmomatic,args=(tumor_fastq_path_first,tumor_fastq_path_second,trimmomatic_path,adapter_path,tumor_fastq_prefix,logfile_out_fold,"tumor",CPU,))
	processes_0.append(q1)
	#for p in processes_0:
	#	p.daemon = True
	#	p.start()
	#for p in processes_0:
	#	p.join()
		
	print "start stage 1"
	processes_1=[]
	if hla_str=="None":
		d1=multiprocessing.Process(target=hlatyping,args=(tumor_fastq_path_first,tumor_fastq_path_second,opitype_fold,opitype_out_fold,opitype_ext,prefix,))
 		#processes_1.append(d1)		
 	else:
 		print "hla type provided!"
	d2=multiprocessing.Process(target=mapping_qc_gatk_preprocess,args=(tumor_fastq_clean_first,tumor_fastq_clean_second,CPU,STAR_INDEX,alignment_out_fold,prefix,REFERENCE,STAR_path,java_picard_path,GATK_path,dbsnp138_path,OneKG_path,mills_path,GTF_path,))
 	#processes_1.append(d2)
 	d3=multiprocessing.Process(target=kallisto_expression,args=(tumor_fastq_clean_first,tumor_fastq_clean_second,kallisto_path,kallisto_out_fold,prefix,kallisto_cdna_path,logfile_out_fold,))
 	processes_1.append(d3)
 	#for p in processes_1:
	#	p.daemon = True
	#	p.start()
	#for p in processes_1:
	#	p.join()
	if hla_str=="None":	
		hla_str=open(opitype_out_fold+'/'+prefix+"_optitype_hla_type").readlines()[0]
	print "start stage 2"
	processes_2=[]
	h1=multiprocessing.Process(target=GATK_hp,args=(GATK_path,REFERENCE,alignment_out_fold,prefix,CPU,dbsnp138_path,somatic_mutation_fold,vcftools_path,vep_path,vep_cache,netmhc_out_fold,pos_1000G_file_path))
 	processes_2.append(h1)
 	#h2=multiprocessing.Process(target=varscan_snv_calling,args=(somatic_mutation_fold,alignment_out_fold,prefix,REFERENCE,vep_cache,samtools_path,varscan_path,vep_path,netmhc_out_fold,))
 	#processes_2.append(h2)
 	#h3=multiprocessing.Process(target=varscan_indel_calling,args=(varscan_indel_fold,alignment_out_fold,prefix,REFERENCE,vep_cache,samtools_path,varscan_path,vep_path,netmhc_out_fold,))
 	#processes_2.append(h3)
 	for p in processes_2:
		p.daemon = True
		p.start()
	for p in processes_2:
		p.join()	
	processes_3=[]
  	l1=multiprocessing.Process(target=varscan_neo,args=(snv_fasta_file,hla_str,snv_netmhc_out_file,netmhc_out_fold,split_num,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_out_fold,netMHCpan_path,))
 	processes_3.append(l1)
   	#l2=multiprocessing.Process(target=indel_neo,args=(prefix,netmhc_out_fold,indel_fasta_file,hla_str,indel_netmhc_out_file,split_num,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_out_fold,netMHCpan_path,))
 	#processes_3.append(l2)		
 	#for p in processes_3:
	#	p.daemon = True
	#	p.start()
	#for p in processes_3:
	#	p.join()

	processes_4=[]
  	r1=multiprocessing.Process(target=immunogenicity_score_calculate,args=(final_neo_file,gmm_classification_file,immunogenicity_score_ranking,immunogenicity_gmm_score_ranking,))
 	processes_4.append(r1)
 	#for p in processes_4:
	#	p.daemon = True
	#	p.start()
	#for p in processes_4:
	#	p.join()






