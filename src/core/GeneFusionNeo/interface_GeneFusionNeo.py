import os
from GeneFusionNeoprocessor import *
import shutil
def GF(opts):
	config_file=opts.Config_file
	f = open(config_file,"r"); # open the configure file
	config_list = {}; # define a hash table
	while 1:
		text = f.readline();
		if text == "":
			break;
		str = text.split(); # split the string
		config_list[str[0]] = str[1]; # assignment the config_list
	f.close()
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
	vep_cache=config_list["vep_cache_path"]
	vep_path=config_list["vep_path"]
	hla_str=config_list["hla_str"]
	split_num=10000
	binding_aff_cutoff=int(config_list["binding_aff_cutoff"])
	fpkm_cutoff=int(config_list["fpkm_cutoff"])
	netMHCpan_path=config_list["netMHCpan_path"]
	gene_fusion_fold_path=output_fold + '/' + "gene_fusion"
	netmhc_out_fold=output_fold + '/' + 'netmhc'
	netctl_out_fold=output_fold + '/' + 'netctl'
	eric_db_path=config_list["eric_db_path"]
	eric_path=config_list["eric_path"]
	genefuison_fasta_file=netmhc_out_fold+'/'+prefix+'_gene_fusion.fasta'
	genefuison_netmhc_out_file=netmhc_out_fold+'/'+prefix+'_genefusion_netmhc.txt'
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

	#####check output directory###
	print "check output directory"
	if not os.path.exists(output_fold):
		os.mkdir(output_fold)	
	if not os.path.exists(netmhc_out_fold):
		os.mkdir(netmhc_out_fold)
	if not os.path.exists(netctl_out_fold):
		os.mkdir(netctl_out_fold)	
	if not os.path.exists(gene_fusion_fold_path):
		os.mkdir(gene_fusion_fold_path)	
	print "start stage 1"
	processes_1=[]
	if hla_str=="None":
		d1=multiprocessing.Process(target=hlatyping,args=(tumor_fastq_path_first,tumor_fastq_path_second,opitype_fold,opitype_out_fold,opitype_ext,prefix,))
 		processes_1.append(d1)
 		hla_str=open(opitype_out_fold+'/'+prefix+"_optitype_hla_type").readlines()[0]		
 	else:
 		print "hla type provided!"
 	gene_fusion(tumor_fastq_path_first,tumor_fastq_path_second,prefix,gene_fusion_fold_path,eric_db_path,eric_path,CPU,netmhc_out_fold,genefuison_fasta_file,hla_str,genefuison_netmhc_out_file,split_num,netMHCpan_path,binding_aff_cutoff,fpkm_cutoff)


