import os
from pairendMDNAprocessor import *
import shutil
import yaml
import time
import warnings
warnings.filterwarnings("ignore")
def PEMD(opts):
	base_dir=os.getcwd()
	config_file=opts.Config_file
	f=open(config_file)
	config_list=yaml.load(f)
	#######read and parse parameter
	print "Read and parse parameters..."
	output_fold_abs=config_list["output_fold"]
	if os.path.isabs(output_fold_abs):
		output_fold=output_fold_abs
	else:
		output_fold=base_dir+'/'+output_fold_abs
	pTuneos_bin_path="bin"
	normal_fastq_path_first=config_list["normal_fastq_path_first"]
	normal_fastq_path_second=config_list["normal_fastq_path_second"]
	tumor_fastq_path_first=config_list["tumor_fastq_path_first"]
	tumor_fastq_path_second=config_list["tumor_fastq_path_second"]
	opitype_fold=config_list["opitype_fold"]
	opitype_out_fold=output_fold + '/' + 'hlatyping'
	opitype_ext=pTuneos_bin_path+'/optitype_ext.py'
	prefix=config_list["sample_name"]
	CPU=config_list["thread_number"]
	vcftools_path="software/vcftools"
	REFERENCE=base_dir + "/" + "database/Fasta/human.fasta"
	BWA_INDEX=base_dir + "/" + "database/Fasta/human.fasta"
	alignment_out_fold=output_fold + '/' + 'alignments'
	bwa_path="software/bwa"
	samtools_path=base_dir + '/' + "software/samtools"
	java_picard_path="software/picard.jar"
	GATK_path="software/GenomeAnalysisTK.jar"
	dbsnp138_path="database/VCF_annotation/dbsnp_138.hg38.vcf.gz"
	OneKG_path="database/VCF_annotation/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
	mills_path="database/VCF_annotation/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
	#cosmic_path="database/VCF_annotation/CosmicCodingMuts_chr_M_sorted.vcf.gz"
	somatic_mutation_fold=output_fold + '/' + 'somatic_mutation'
	vep_cache=config_list["vep_cache_path"]
	varscan_path=base_dir + '/' + "software/VarScan.v2.4.2.jar"
	vep_path=config_list["vep_path"]
	clean_fastq_fold=output_fold + '/' + 'clean_fastq'
	netmhc_out_fold=output_fold + '/' + 'netmhc'
	strelka_out_fold=output_fold + '/' + 'strelka'
	logfile_out_fold=output_fold + '/' + 'logfile'
	bamstat_out_fold=output_fold + '/' + 'bamstat'
	indel_fasta_file=netmhc_out_fold+'/'+prefix+'_indel.fasta'
	hla_str=config_list["hla_str"]
	indel_netmhc_out_file=netmhc_out_fold+'/'+prefix+'_indel_netmhc.tsv'
	split_num=200
	human_peptide_path="database/Protein/human.pep.all.fa"
	binding_fc_aff_cutoff=int(config_list["binding_fc_aff_cutoff"])
	binding_aff_cutoff=int(config_list["binding_aff_cutoff"])
	fpkm_cutoff=int(config_list["fpkm_cutoff"])
	netctl_out_fold=output_fold + '/' + 'netctl'
	netMHCpan_path=config_list["netMHCpan_path"]
	varscan_copynumber_fold=output_fold + '/' + 'copynumber_profile'
	snv_fasta_file=netmhc_out_fold+'/'+prefix+'_snv.fasta'
	snv_netmhc_out_file=netmhc_out_fold+'/'+prefix+'_snv_netmhc.tsv'
	pyclone_fold=output_fold + '/' + 'pyclone'
	pyclone_path=config_list["pyclone_path"]
	strelka_path=config_list["strelka_path"]
	rna_fastq_1_path=config_list["tumor_rna_fastq_1"]
	rna_fastq_2_path=config_list["tumor_rna_fastq_2"]
	kallisto_path="software/kallisto"
	tumor_depth_cutoff=config_list["tumor_depth_cutoff"]
	tumor_vaf_cutoff=config_list["tumor_vaf_cutoff"]
	normal_vaf_cutoff=config_list["normal_vaf_cutoff"]
	kallisto_out_fold=output_fold + '/' + 'expression'
	kallisto_cdna_path="database/Protein/human.cdna.all.fa"
	snv_final_neo_file=netctl_out_fold + '/' + prefix + '_pyclone_neo.tsv'
	indel_final_neo_file=netctl_out_fold + '/' + prefix + '_indel_netctl_concact.tsv'
	trimmomatic_path="software/trimmomatic-0.36.jar"
	netchop_path="software/netchop"
	adapter_path="software/TruSeq3-PE.fa"
	tumor_fastq_prefix=clean_fastq_fold + '/' + prefix + "_tumor_clean.fq.gz"
	normal_fastq_prefix=clean_fastq_fold + '/' + prefix + "_normal_clean.fq.gz"
	tumor_fastq_dir=clean_fastq_fold
	tumor_fastq_clean_first=clean_fastq_fold + '/' + prefix + "_tumor_clean_1P.fq.gz"
	tumor_fastq_clean_second=clean_fastq_fold + '/' + prefix + "_tumor_clean_2P.fq.gz"
	normal_fastq_dir=clean_fastq_fold
	normal_fastq_clean_first=clean_fastq_fold + '/' + prefix + "_normal_clean_1P.fq.gz"
	normal_fastq_clean_second=clean_fastq_fold + '/' + prefix + "_normal_clean_2P.fq.gz"
	iedb_file="train_model/iedb.fasta"
	peptide_length=config_list["peptide_length"]
	cf_hy_model_9="train_model/cf_hy_9_model.m"
	cf_hy_model_10="train_model/cf_hy_10_model.m"
	cf_hy_model_11="train_model/cf_hy_11_model.m"
	RF_model="train_model/RF_train_model.m"
	driver_gene_path="software/DriveGene.tsv"
	snv_neo_model_file=netctl_out_fold + '/' + prefix + '_snv_neo_model.tsv'
	snv_blastp_tmp_file=netctl_out_fold + '/' + prefix + '_snv_blastp_tmp.tsv'
	snv_blastp_out_tmp_file=netctl_out_fold + '/' + prefix + '_snv_blastp_out_tmp.tsv'
	snv_netMHCpan_pep_tmp_file=netctl_out_fold + '/' + prefix + '_snv_netMHCpan_pep_tmp.tsv'
	snv_netMHCpan_ml_out_tmp_file=netctl_out_fold + '/' + prefix + '_snv_netMHCpan_ml_out_tmp.tsv'
	indel_neo_model_file=netctl_out_fold + '/' + prefix + '_indel_neo_model.tsv'
	indel_blastp_tmp_file=netctl_out_fold + '/' + prefix + '_indel_blastp_tmp.tsv'
	indel_blastp_out_tmp_file=netctl_out_fold + '/' + prefix + '_indel_blastp_out_tmp.tsv'
	indel_netMHCpan_pep_tmp_file=netctl_out_fold + '/' + prefix + '_indel_netMHCpan_pep_tmp.tsv'
	indel_netMHCpan_ml_out_tmp_file=netctl_out_fold + '/' + prefix + '_indel_netMHCpan_ml_out_tmp.tsv'
	blast_db_path="database/Protein/peptide_database/peptide"
	keep_tmp=int(config_list["tmp"])
	fragment_length=config_list["fragment_length"]
	fragment_SD=config_list["fragment_SD"]
	time.sleep(5)
	print "Read and parse parameters...  OK"
	print "Check reference file path and input file path..."
	if os.path.exists(rna_fastq_1_path) and os.path.exists(rna_fastq_2_path):
		exp_file=output_fold + '/' + "expression/abundance.tsv"
	else:
		exp_file=config_list["expression_file"]
	#####check input file,tool path and reference file#####
	if os.path.exists(normal_fastq_path_first) and os.path.exists(normal_fastq_path_second) and os.path.exists(tumor_fastq_path_first) and os.path.exists(tumor_fastq_path_second) and os.path.exists(rna_fastq_1_path):
		print "Check all fastq file...  OK"
	else:
		print "Please check your input fastq file!"
		os._exit(1)
	if os.path.exists(opitype_fold):
		print "Check opitype path...  OK"
	else:
		print "Please check your opitype path!"
		os._exit(1)
	if os.path.exists(bwa_path):
		print "Check bwa path...  OK"
	else:
		print "Please check your bwa path!"
		os._exit(1)	
	if os.path.exists(kallisto_path):
		print "Check kallisto path...  OK"
	else:
		print "Please check your kallisto path!"
		os._exit(1)	
	if os.path.exists(samtools_path):
		print "Check samtools path...  OK"
	else:
		print "Please check your samtools path!"
		os._exit(1)	
	if os.path.exists(java_picard_path):
		print "Check picard path...  OK"
	else:
		print "Please check your picard path!"
		os._exit(1)	
	if os.path.exists(GATK_path):
		print "Check GATK path...  OK"
	else:
		print "Please check your GATK path!"
		os._exit(1)	
	if os.path.exists(dbsnp138_path):
		print "Check dbsnp138 file path...  OK"
	else:
		print "Please check your dbsnp138 file path!"
		os._exit(1)	
	if os.path.exists(OneKG_path):
		print "Check OneKG file path...  OK"
	else:
		print "Please check your OneKG file path!"
		os._exit(1)	
	if os.path.exists(mills_path):
		print "Check mills file path...  OK"
	else:
		print "Please check your mills file path!"
		os._exit(1)	
	if os.path.exists(vep_path):
		print "Check vep path...  OK"
	else:
		print "Please check your vep path!"
		os._exit(1)	
	if os.path.exists(vep_cache):
		print "Check vep cache path... OK"
	else:
		print "Please check your vep cache path!"
		os._exit(1)		
	if os.path.exists(REFERENCE):
		print "Check REFERENCE file path...  OK"
	else:
		print "Please check your REFERENCE file path!"
		os._exit(1)	
	if os.path.exists(kallisto_cdna_path): 
		print "Kallisto cdna reference path:%s"%(kallisto_cdna_path)
	else:
		print "ERROR: no kallisto cdna reference file"
		os._exit(1)
	time.sleep(5)
	#####check output directory###
	print "Check output directory"
	if not os.path.exists(output_fold):
		os.mkdir(output_fold)
	if not os.path.exists(somatic_mutation_fold):
		os.mkdir(somatic_mutation_fold)
	if not os.path.exists(netmhc_out_fold):
		os.mkdir(netmhc_out_fold)
	if not os.path.exists(netctl_out_fold):
		os.mkdir(netctl_out_fold)
	if not os.path.exists(varscan_copynumber_fold):
		os.mkdir(varscan_copynumber_fold)
	if not os.path.exists(pyclone_fold):
		os.mkdir(pyclone_fold)
	if not os.path.exists(alignment_out_fold):
		os.mkdir(alignment_out_fold)
	if not os.path.exists(kallisto_out_fold):
		os.mkdir(kallisto_out_fold)
	if not os.path.exists(logfile_out_fold):
		os.mkdir(logfile_out_fold)	
	if not os.path.exists(bamstat_out_fold):
		os.mkdir(bamstat_out_fold)
	if not os.path.exists(clean_fastq_fold):
		os.mkdir(clean_fastq_fold)
	if not os.path.exists(opitype_out_fold):
		os.mkdir(opitype_out_fold)	
	print "Check output directory paths...  OK"	
	time.sleep(5)
	print "Start stage 0: sequence quality control."
	processes_0=[]
	q1=multiprocessing.Process(target=read_trimmomatic,args=(tumor_fastq_path_first,tumor_fastq_path_second,trimmomatic_path,adapter_path,tumor_fastq_prefix,logfile_out_fold,"tumor",CPU,))
	processes_0.append(q1)
	q2=multiprocessing.Process(target=read_trimmomatic,args=(normal_fastq_path_first,normal_fastq_path_second,trimmomatic_path,adapter_path,normal_fastq_prefix,logfile_out_fold,"normal",CPU,))
	processes_0.append(q2)
	#for p in processes_0:
	#	p.daemon = True
	#	p.start()
	#for p in processes_0:
	#	p.join()
	print "Stage 0 finished!"
	print "Start stage 1: hla typing, sequence mapping and expression profiling!"
	processes_1=[]
	if hla_str=="None":
		d1=multiprocessing.Process(target=hlatyping,args=(tumor_fastq_path_first,tumor_fastq_path_second,opitype_fold,opitype_out_fold,opitype_ext,prefix,logfile_out_fold,))
 		#processes_1.append(d1)
 	else:
 		print "hla type provided!"
 	d2=multiprocessing.Process(target=mapping_qc_gatk_preprocess,args=(normal_fastq_clean_first,normal_fastq_clean_second,'normal',CPU,BWA_INDEX,alignment_out_fold,prefix,REFERENCE,bwa_path,samtools_path,java_picard_path,GATK_path,dbsnp138_path,OneKG_path,mills_path,logfile_out_fold,bamstat_out_fold,))
 	#processes_1.append(d2)
 	d3=multiprocessing.Process(target=mapping_qc_gatk_preprocess,args=(tumor_fastq_path_first,tumor_fastq_path_second,'tumor',CPU,BWA_INDEX,alignment_out_fold,prefix,REFERENCE,bwa_path,samtools_path,java_picard_path,GATK_path,dbsnp138_path,OneKG_path,mills_path,logfile_out_fold,bamstat_out_fold,))
 	#processes_1.append(d3)
 	if os.path.exists(rna_fastq_1_path):
 		d4=multiprocessing.Process(target=kallisto_expression,args=(rna_fastq_1_path,rna_fastq_2_path,kallisto_path,kallisto_out_fold,prefix,kallisto_cdna_path,logfile_out_fold,CPU,fragment_length,fragment_SD,))
 		processes_1.append(d4)
 	else:
 		print "RNA sequence not found, the kallisto will not be run!"
 	for p in processes_1:
		p.daemon = True
		p.start()
	for p in processes_1:
		p.join()
	print 'Stage 1 finished!'
	print 'Start stage 2: mutation calling using Mutect2, strelka and varscan; copynumber profiling.'
	processes_2=[]
	h0=multiprocessing.Process(target=GATK_mutect2,args=(GATK_path,REFERENCE,alignment_out_fold,prefix,CPU,dbsnp138_path,somatic_mutation_fold,vcftools_path,vep_path,vep_cache,netmhc_out_fold,tumor_depth_cutoff,tumor_vaf_cutoff,normal_vaf_cutoff,pTuneos_bin_path,human_peptide_path,logfile_out_fold))
	processes_2.append(h0)
	h1=multiprocessing.Process(target=varscan_somatic_caling_drift,args=(somatic_mutation_fold,alignment_out_fold,prefix,REFERENCE,vep_cache,samtools_path,varscan_path,vep_path,netmhc_out_fold,logfile_out_fold,))
	processes_2.append(h1)
	h2=multiprocessing.Process(target=indel_calling_drift,args=(strelka_out_fold,strelka_path,alignment_out_fold,prefix,REFERENCE,vep_cache,netmhc_out_fold,CPU,vep_path,pTuneos_bin_path,logfile_out_fold,human_peptide_path,))
	processes_2.append(h2)
	h3=multiprocessing.Process(target=varscan_copynumber_calling,args=(varscan_copynumber_fold,prefix,alignment_out_fold,REFERENCE,samtools_path,varscan_path,logfile_out_fold))
	processes_2.append(h3)
	#for p in processes_2:
	#	p.daemon = True
	#	p.start()
	#for p in processes_2:
	#	p.join()
	#if hla_str=="None":
 	#	hla_str=open(opitype_out_fold+'/'+prefix+"_optitype_hla_type").readlines()[0]
	print 'Stage 2 finished!'
	print 'Start stage 3: neoantigens identification.'
	processes_3=[]
	t2=multiprocessing.Process(target=varscan_neo,args=(snv_fasta_file,hla_str,driver_gene_path,snv_netmhc_out_file,netmhc_out_fold,split_num,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_out_fold,netMHCpan_path,pTuneos_bin_path,peptide_length,netchop_path,))
	processes_3.append(t2)
	t3=multiprocessing.Process(target=indel_neo,args=(somatic_mutation_fold,prefix,vep_cache,netmhc_out_fold,vep_path,indel_fasta_file,hla_str,driver_gene_path,indel_netmhc_out_file,split_num,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_out_fold,netMHCpan_path,pTuneos_bin_path,peptide_length,netchop_path,REFERENCE,human_peptide_path,))
	processes_3.append(t3)
	#for p in processes_3:
	#	p.daemon = True
	#	p.start()
	#for p in processes_3:
	#	p.join()
	print "Stage 3 finished."
	print 'Start stage 4: mutation clonal cellularity calculation.'
	processes_4=[]
	l1=multiprocessing.Process(target=pyclone_annotation,args=(somatic_mutation_fold,varscan_copynumber_fold,prefix,pyclone_fold,netctl_out_fold,pyclone_path,logfile_out_fold,pTuneos_bin_path,))
	processes_4.append(l1)	
	#for p in processes_4:
	#	p.daemon = True
	#	p.start()
	#for p in processes_4:
	#	p.join()
	print 'Stage 4 finished.'
	print 'Start stage 5: neoantigen filtering using Pre&RecNeo model and refined immunogenicity score scheme.'
	processes_5=[]
	r1=multiprocessing.Process(target=InVivoModelAndScoreSNV,args=(snv_final_neo_file,cf_hy_model_9,cf_hy_model_10,cf_hy_model_11,RF_model,snv_neo_model_file,snv_blastp_tmp_file,snv_blastp_out_tmp_file,snv_netMHCpan_pep_tmp_file,snv_netMHCpan_ml_out_tmp_file,iedb_file,blast_db_path,))
	processes_5.append(r1)
	r2=multiprocessing.Process(target=InVivoModelAndScoreINDEL,args=(indel_final_neo_file,cf_hy_model_9,cf_hy_model_10,cf_hy_model_11,RF_model,indel_neo_model_file,indel_blastp_tmp_file,indel_blastp_out_tmp_file,indel_netMHCpan_pep_tmp_file,indel_netMHCpan_ml_out_tmp_file,iedb_file,blast_db_path,))
	processes_5.append(r2)
	#for p in processes_5:
	#	p.daemon = True
	#	p.start()
	#for p in processes_5:
	#	p.join()		
	print 'Stage 5 finished.'
'''
	if keep_tmp==0:
		if os.path.exists(snv_blastp_tmp_file):
			os.remove(snv_blastp_tmp_file)
		if os.path.exists(snv_blastp_out_tmp_file):
			os.remove(snv_blastp_out_tmp_file)
		if os.path.exists(snv_netMHCpan_pep_tmp_file):
			os.remove(snv_netMHCpan_pep_tmp_file)
		if os.path.exists(snv_netMHCpan_ml_out_tmp_file):
			os.remove(snv_netMHCpan_ml_out_tmp_file)
		if os.path.exists(indel_blastp_tmp_file):
			os.remove(indel_blastp_tmp_file)
		if os.path.exists(indel_blastp_out_tmp_file):
			os.remove(indel_blastp_out_tmp_file)
		if os.path.exists(indel_netMHCpan_pep_tmp_file):
			os.remove(indel_netMHCpan_pep_tmp_file)
		if os.path.exists(indel_netMHCpan_ml_out_tmp_file):
			os.remove(indel_netMHCpan_ml_out_tmp_file)
		if os.path.exists(snv_neo_model_file):
			if os.path.getsize(snv_neo_model_file):
				shutil.rmtree(alignment_out_fold)
				shutil.rmtree(clean_fastq_fold)
				shutil.rmtree(somatic_mutation_fold)
				shutil.rmtree(varscan_copynumber_fold)
				shutil.rmtree(netmhc_out_fold)
				shutil.rmtree(pyclone_fold)
				shutil.rmtree(kallisto_out_fold)
				shutil.rmtree(opitype_out_fold)
				shutil.rmtree(strelka_out_fold)
				shutil.rmtree(logfile_out_fold)
	else:
		print "Keep all tmporal files!"
	print "ALL finished! Please check result files 'snv_neo_model.tsv' and 'indel_neo_model.tsv' in netctl fold"
'''