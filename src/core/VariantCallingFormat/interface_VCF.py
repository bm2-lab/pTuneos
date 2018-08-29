import os,multiprocessing
from VCFprocessor import *
import shutil
import yaml
def Vcf(opts):
	config_file=opts.Config_file
	f=open(config_file)
	config_list=yaml.load(f)
	#######read and parse parameter
	print "read and parse parameter."
	output_fold=config_list["output_fold"]
	itunes_bin_path=config_list["itunes_bin_path"]
	os.system("export iTuNES_BIN_PATH=%s"%itunes_bin_path)
	vcf_file=config_list["vcf_file"]
	somatic_out_fold=output_fold + '/' + 'somatic_mutation'
	prefix=config_list["sample_name"]
	REFERENCE=config_list["reference_path"]
	tumor_depth_cutoff=config_list["tumor_depth_cutoff"]
	tumor_vaf_cutoff=config_list["tumor_vaf_cutoff"]
	normal_vaf_cutoff=config_list["normal_vaf_cutoff"]
	vep_cache=config_list["vep_cache_path"]
	vep_path=config_list["vep_path"]
	netmhc_out_fold=output_fold + '/' + 'netmhc'
	indel_fasta_file=netmhc_out_fold+'/'+prefix+'_indel.fasta'
	hla_str=config_list["hla_str"]
	split_num=200
	exp_file=config_list["expression_file"]
	binding_fc_aff_cutoff=int(config_list["binding_fc_aff_cutoff"])
	binding_aff_cutoff=int(config_list["binding_aff_cutoff"])
	fpkm_cutoff=int(config_list["fpkm_cutoff"])
	netctl_out_fold=output_fold + '/' + 'netctl'
	netMHCpan_path=config_list["netMHCpan_path"]
	snv_fasta_file=netmhc_out_fold+'/'+prefix+'_snv.fasta'
	snv_netmhc_out_file=netmhc_out_fold+'/'+prefix+'_snv_netmhc.txt'
	indel_netmhc_out_file=netmhc_out_fold+'/'+prefix+'_indel_netmhc.txt'
	vcftools_path=config_list["vcftools_path"]
	peptide_length=config_list["peptide_length"]
	coverage=int(config_list["coverage_cutoff"])
	pyclone_fold=output_fold + '/' + 'pyclone'
	pyclone_path=config_list["pyclone_path"]
	copynumber_profile=config_list["copynumber_profile"]
	tumor_cellularity=float(config_list["tumor_cellularity"])
	snv_final_neo_file=netctl_out_fold + '/' + prefix + '_pyclone_neo.txt'
	snv_gmm_classification_file=netctl_out_fold + '/' + prefix + '_snv_gmm.png'
	snv_immunogenicity_gmm_all_score_ranking=netctl_out_fold + '/' + prefix + '_snv_gmm_all_score.txt'
	snv_immunogenicity_gmm_pos_score_ranking=netctl_out_fold + '/' + prefix + '_snv_gmm_pos_score.txt'
	snv_immunogenicity_bioactive_score_ranking=netctl_out_fold + '/' + prefix + '_snv_bioactive_score.txt'
	indel_gmm_classification_file=netctl_out_fold + '/' + prefix + '_indel_gmm.png'
	indel_immunogenicity_gmm_all_score_ranking=netctl_out_fold + '/' + prefix + '_indel_gmm_all_score.txt'
	indel_immunogenicity_gmm_pos_score_ranking=netctl_out_fold + '/' + prefix + '_indel_gmm_pos_score.txt'
	indel_immunogenicity_bioactive_score_ranking=netctl_out_fold + '/' + prefix + '_indel_bioactive_score.txt'
	iedb_file=config_list["iedb_file"]
	mhc_pos_file=config_list["mhc_pos_file"]
	mhc_neg_file=config_list["mhc_neg_file"]
	model_train_file=config_list['model_train_file']
	snv_neo_model_file=netctl_out_fold + '/' + prefix + '_snv_neo_model.txt'
	snv_blastp_tmp_file=netctl_out_fold + '/' + prefix + '_snv_blastp_tmp.txt'
	snv_blastp_out_tmp_file=netctl_out_fold + '/' + prefix + '_snv_blastp_out_tmp.txt'
	snv_netMHCpan_pep_tmp_file=netctl_out_fold + '/' + prefix + '_snv_netMHCpan_pep_tmp.txt'
	snv_netMHCpan_ml_out_tmp_file=netctl_out_fold + '/' + prefix + '_snv_netMHCpan_ml_out_tmp.txt'
	indel_neo_model_file=netctl_out_fold + '/' + prefix + '_indel_neo_model.txt'
	indel_blastp_tmp_file=netctl_out_fold + '/' + prefix + '_indel_blastp_tmp.txt'
	indel_blastp_out_tmp_file=netctl_out_fold + '/' + prefix + '_indel_blastp_out_tmp.txt'
	indel_netMHCpan_pep_tmp_file=netctl_out_fold + '/' + prefix + '_indel_netMHCpan_pep_tmp.txt'
	indel_netMHCpan_ml_out_tmp_file=netctl_out_fold + '/' + prefix + '_indel_netMHCpan_ml_out_tmp.txt'
	blast_db_path=config_list['blast_db_path']
	print vcftools_path
	#####check input file,tool path and reference file#####
	if os.path.exists(vcf_file):
		print "check vcf file done."
	else:
		print "please check your input vcf file!"
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
	if exp_file!="no_exp" and os.path.exists(exp_file):
		print "check expression file done."
	elif exp_file=="no_exp":
		print "no expression file provided."
	else:
		print "please check your expression file path!"
		os._exit(1)	
	if os.path.exists(REFERENCE):
		print "check REFERENCE file path done."
	else:
		print "please check your REFERENCE file path!"
		os._exit(1)	

	#####check output directory###
	print "check output directory"
	if not os.path.exists(output_fold):
		os.mkdir(output_fold)
	if not os.path.exists(somatic_out_fold):
		os.mkdir(somatic_out_fold)
	if not os.path.exists(netmhc_out_fold):
		os.mkdir(netmhc_out_fold)
	if not os.path.exists(netctl_out_fold):
		os.mkdir(netctl_out_fold)
	if hla_str=="None":
		print "please provied hla type, seperate by comma."		
 	else:
 		print "hla type provided!"
  	processes_0=[]
	h1=multiprocessing.Process(target=VCF_process,args=(prefix,vcf_file,somatic_out_fold,vcftools_path,vep_path,vep_cache,netmhc_out_fold,tumor_depth_cutoff,tumor_vaf_cutoff,normal_vaf_cutoff,))
 	processes_0.append(h1)
 	for p in processes_0:
		p.daemon = True
		p.start()
	for p in processes_0:
		p.join()
	print "start stage 1"
 	processes_1=[]
	d1=multiprocessing.Process(target=snv_neo,args=(snv_fasta_file,hla_str,snv_netmhc_out_file,netmhc_out_fold,split_num,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_out_fold,netMHCpan_path,peptide_length,))
 	processes_1.append(d1)
 	d2=multiprocessing.Process(target=indel_neo,args=(indel_fasta_file,somatic_out_fold,hla_str,indel_netmhc_out_file,split_num,netMHCpan_path,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_out_fold,netmhc_out_fold,peptide_length,))
 	processes_1.append(d2)
 	#q.put('tumor_qc')
 	for p in processes_1:
		p.daemon = True
		p.start()
	for p in processes_1:
		p.join()
	print "start stage 2"
 	processes_2=[]
	m1=multiprocessing.Process(target=pyclone_annotation,args=(somatic_out_fold,copynumber_profile,tumor_cellularity,prefix,pyclone_fold,netctl_out_fold,coverage,pyclone_path,))
 	processes_2.append(m1)
 	for p in processes_2:
		p.daemon = True
		p.start()
	for p in processes_2:
		p.join()
	print "start stage 2"
	processes_3=[]
	r1=multiprocessing.Process(target=InVivoModelAndScoreSNV,args=(mhc_pos_file,mhc_neg_file,snv_final_neo_file,model_train_file,snv_neo_model_file,snv_blastp_tmp_file,snv_blastp_out_tmp_file,snv_netMHCpan_pep_tmp_file,snv_netMHCpan_ml_out_tmp_file,iedb_file,blast_db_path,snv_immunogenicity_gmm_all_score_ranking,snv_immunogenicity_gmm_pos_score_ranking,snv_gmm_classification_file,snv_immunogenicity_bioactive_score_ranking,))
	processes_3.append(r1)
	#r2=multiprocessing.Process(target=InVivoModelAndScoreINDEL,args=(mhc_pos_file,mhc_neg_file,indel_final_neo_file,model_train_file,indel_neo_model_file,indel_blastp_tmp_file,indel_blastp_out_tmp_file,indel_netMHCpan_pep_tmp_file,indel_netMHCpan_ml_out_tmp_file,iedb_file,blast_db_path,indel_immunogenicity_gmm_all_score_ranking,indel_immunogenicity_gmm_pos_score_ranking,indel_gmm_classification_file,indel_immunogenicity_bioactive_score_ranking,))
	#processes_3.append(r2)
	for p in processes_3:
		p.daemon = True
		p.start()
	for p in processes_3:
		p.join()		
