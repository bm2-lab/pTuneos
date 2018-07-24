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
	prefix=config_list["sample_name"]
	REFERENCE=config_list["reference_path"]
	vcf_fold=output_fold + '/' + 'vcf'
	vep_cache=config_list["vep_cache_path"]
	vep_path=config_list["vep_path"]
	netmhc_out_fold=output_fold + '/' + 'netmhc'
	indel_fasta_file=netmhc_out_fold+'/'+prefix+'_indel.fasta'
	hla_str=config_list["hla_str"]
	indel_netmhc_out_file=netmhc_out_fold+'/'+prefix+'_indel_netmhc.txt'
	split_num=200
	exp_file=config_list["expression_file"]
	binding_fc_aff_cutoff=int(config_list["binding_fc_aff_cutoff"])
	binding_aff_cutoff=int(config_list["binding_aff_cutoff"])
	fpkm_cutoff=int(config_list["fpkm_cutoff"])
	netctl_out_fold=output_fold + '/' + 'netctl'
	netMHCpan_path=config_list["netMHCpan_path"]
	snv_fasta_file=netmhc_out_fold+'/'+prefix+'_snv.fasta'
	snv_netmhc_out_file=netmhc_out_fold+'/'+prefix+'_snv_netmhc.txt'
	reference_build=config_list["reference_build"]
	crossmap_path=config_list["crossmap_path"]
	vcftools_path=config_list["vcftools_path"]
	print vcftools_path
	snv_vcf_file=vcf_fold + '/' + prefix + '_SNVs_only.recode.vcf'
	indel_vcf_file=vcf_fold + '/' + prefix + '_INDELs_only.recode.vcf'
	liftover_file=config_list["liftover_file"]
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
	if not os.path.exists(vcf_fold):
		os.mkdir(vcf_fold)
	if not os.path.exists(netmhc_out_fold):
		os.mkdir(netmhc_out_fold)
	if not os.path.exists(netctl_out_fold):
		os.mkdir(netctl_out_fold)
	if hla_str=="None":
		print "please provied hla type, seperate by comma."		
 	else:
 		print "hla type provided!"
 	#if reference_build=="GRch37":
 	#	hg_liftover(crossmap_path,liftover_file,vcf_file,REFERENCE,vcf_fold,vcftools_path,prefix)
 	#else:
 	#	hg_no_liftover(vcf_file,vcf_fold,vcftools_path,prefix)
 	processes_1=[]
	d1=multiprocessing.Process(target=snv_neo,args=(vep_path,snv_vcf_file,vep_cache,vcf_fold,netmhc_out_fold,snv_fasta_file,snv_netmhc_out_file,split_num,netMHCpan_path,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,hla_str,netctl_out_fold,))
 	processes_1.append(d1)
 	#d2=multiprocessing.Process(target=indel_neo,args=(vep_path,indel_vcf_file,vep_cache,vcf_fold,netmhc_out_fold,indel_fasta_file,hla_str,indel_netmhc_out_file,split_num,netMHCpan_path,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_out_fold,))
 	#processes_1.append(d2)
 	#q.put('tumor_qc')
 	for p in processes_1:
		p.daemon = True
		p.start()
	for p in processes_1:
		p.join()
	print "start stage 2"

	
