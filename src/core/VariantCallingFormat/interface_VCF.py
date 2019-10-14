import os,multiprocessing
from VCFprocessor import *
import shutil
import yaml
import time
import warnings
warnings.filterwarnings("ignore")
def Vcf(opts):
	base_dir=os.getcwd()
	config_file=opts.Config_file
	f=open(config_file)
	config_list=yaml.load(f)
	#######read and parse parameter
	print "Start reading and parsing parameter..."
	time.sleep(5)
	output_fold=config_list["output_fold"]
	ptuneos_bin_path="bin"
	vcf_file=config_list["vcf_file"]
	REFERENCE=base_dir + "/" + "database/Fasta/human.fasta"
	somatic_out_fold=output_fold + '/' + 'somatic_mutation'
	logfile_out_fold=output_fold + '/' + 'logfile'
	prefix=config_list["sample_name"]
	tumor_depth_cutoff=config_list["tumor_depth_cutoff"]
	tumor_vaf_cutoff=config_list["tumor_vaf_cutoff"]
	normal_vaf_cutoff=config_list["normal_vaf_cutoff"]
	vep_cache=config_list["vep_cache_path"]
	vep_path=config_list["vep_path"]
	netmhc_out_fold=output_fold + '/' + 'netmhc'
	indel_fasta_file=netmhc_out_fold+'/'+prefix+'_indel.fasta'
	hla_str=config_list["hla_str"]
	split_num=200
	netchop_path="software/netchop"
	human_peptide_path="database/Protein/human.pep.all.fa"
	exp_file=config_list["expression_file"]
	binding_fc_aff_cutoff=int(config_list["binding_fc_aff_cutoff"])
	binding_aff_cutoff=int(config_list["binding_aff_cutoff"])
	fpkm_cutoff=int(config_list["fpkm_cutoff"])
	netctl_out_fold=output_fold + '/' + 'netctl'
	netMHCpan_path=config_list["netMHCpan_path"]
	snv_fasta_file=netmhc_out_fold+'/'+prefix+'_snv.fasta'
	snv_netmhc_out_file=netmhc_out_fold+'/'+prefix+'_snv_netmhc.tsv'
	indel_netmhc_out_file=netmhc_out_fold+'/'+prefix+'_indel_netmhc.tsv'
	vcftools_path="software/vcftools"
	peptide_length=config_list["peptide_length"]
	pyclone_fold=output_fold + '/' + 'pyclone'
	pyclone_path=config_list["pyclone_path"]
	copynumber_profile=config_list["copynumber_profile"]
	tumor_cellularity=float(config_list["tumor_cellularity"])
	snv_final_neo_file=netctl_out_fold + '/' + prefix + '_pyclone_neo.tsv'
	indel_final_neo_file=netctl_out_fold + '/' + prefix + '_indel_netctl_concact.tsv'
	iedb_file="train_model/iedb.fasta"
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
	#####check input file,tool path and reference file#####
	if os.path.exists(vcf_file):
		print "Check inuput mutation vcf file...  OK"
	else:
		print "Please check your input vcf file!"
		os._exit(1)
	if os.path.exists(vep_path):
		print "Check vep path...  OK"
	else:
		print "Please check your vep path!"
		os._exit(1)	
	if os.path.exists(vep_cache):
		print "Check vep cache path...  OK"
	else:
		print "Please check your vep cache path!"
		os._exit(1)	
	if os.path.exists(exp_file):
		print "Check expression file...  OK"
	else:
		print "Please check your expression file path!"
		os._exit(1)		
	time.sleep(5)
	#####check output directory###
	print "Check output directory"
	if not os.path.exists(output_fold):
		os.mkdir(output_fold)
	if not os.path.exists(somatic_out_fold):
		os.mkdir(somatic_out_fold)
	if not os.path.exists(netmhc_out_fold):
		os.mkdir(netmhc_out_fold)
	if not os.path.exists(netctl_out_fold):
		os.mkdir(netctl_out_fold)
	if not os.path.exists(logfile_out_fold):
		os.mkdir(logfile_out_fold)
	if not os.path.exists(pyclone_fold):
		os.mkdir(pyclone_fold)
	if hla_str=="None":
		print "please provied hla type, seperate by comma,eg:HLiA-A02:01,HLA-A01:01,HLA-B15:17,HLA-B13:02,HLA-C07:01,HLA-C06:02"		
 	else:
 		print "Check hla alleles...  OK"
 	print "Start preprocessing VCF file..."
  	processes_0=[]
	h1=multiprocessing.Process(target=VCF_process,args=(prefix,vcf_file,somatic_out_fold,vcftools_path,vep_path,vep_cache,netmhc_out_fold,tumor_depth_cutoff,tumor_vaf_cutoff,normal_vaf_cutoff,ptuneos_bin_path,human_peptide_path,logfile_out_fold,))
 	processes_0.append(h1)
 	for p in processes_0:
		p.daemon = True
		p.start()
	for p in processes_0:
		p.join()
	print "Preprocessing VCF file done!"
	print "Start neoantigen prediction..."
 	processes_1=[]
	d1=multiprocessing.Process(target=snv_neo,args=(snv_fasta_file,hla_str,driver_gene_path,snv_netmhc_out_file,netmhc_out_fold,split_num,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_out_fold,netMHCpan_path,peptide_length,ptuneos_bin_path,netchop_path,))
 	processes_1.append(d1)
 	d2=multiprocessing.Process(target=indel_neo,args=(indel_fasta_file,somatic_out_fold,hla_str,driver_gene_path,indel_netmhc_out_file,split_num,netMHCpan_path,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_out_fold,netmhc_out_fold,peptide_length,ptuneos_bin_path,netchop_path,REFERENCE,human_peptide_path,))
 	processes_1.append(d2)
 	for p in processes_1:
		p.daemon = True
		p.start()
	for p in processes_1:
		p.join()
	print "Neoantigen prediciton done!"

	print "Neoantigen annotation..."
 	processes_2=[]
	m1=multiprocessing.Process(target=pyclone_annotation,args=(somatic_out_fold,copynumber_profile,tumor_cellularity,prefix,pyclone_fold,netctl_out_fold,pyclone_path,ptuneos_bin_path,logfile_out_fold,))
 	processes_2.append(m1)
 	for p in processes_2:
		p.daemon = True
		p.start()
	for p in processes_2:
		p.join()
	print "Neoantigen annotation done!"
	print "Neoantigen filtering using Pre&RecNeo model and refined immunogenicity score scheme."
	processes_3=[]
	r1=multiprocessing.Process(target=InVivoModelAndScoreSNV,args=(snv_final_neo_file,cf_hy_model_9,cf_hy_model_10,cf_hy_model_11,RF_model,snv_neo_model_file,snv_blastp_tmp_file,snv_blastp_out_tmp_file,snv_netMHCpan_pep_tmp_file,snv_netMHCpan_ml_out_tmp_file,iedb_file,blast_db_path,))
	processes_3.append(r1)
	if os.path.exists(indel_final_neo_file):
		r2=multiprocessing.Process(target=InVivoModelAndScoreINDEL,args=(indel_final_neo_file,cf_hy_model_9,cf_hy_model_10,cf_hy_model_11,RF_model,indel_neo_model_file,indel_blastp_tmp_file,indel_blastp_out_tmp_file,indel_netMHCpan_pep_tmp_file,indel_netMHCpan_ml_out_tmp_file,iedb_file,blast_db_path,))
		processes_3.append(r2)
	else:
		print "No neoantigen from Indels is identified!"
	for p in processes_3:
		p.daemon = True
		p.start()
	for p in processes_3:
		p.join()		
	print "All Finished! please check result files 'snv_neo_model.tsv' and 'indel_neo_model.tsv' in netctl fold"
