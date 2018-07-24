import os,sys,time
import multiprocessing
import shutil 
import subprocess
import pandas as pd
import math
from pyper import *
import numpy as np
import math
from sklearn import preprocessing
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from sklearn.semi_supervised import label_propagation
import itertools
from scipy import linalg
import matplotlib as mpl
from sklearn import mixture

def read_trimmomatic(raw_fastq_path_first,raw_fastq_path_second,trimmomatic_path,adapter_path,fastq_prefix,logfile_fold,fastq_type,CPU):
	cmd_trimmomatic="java -jar " + trimmomatic_path + " PE -threads " + str(CPU) + " -phred33 " + raw_fastq_path_first + ' ' + raw_fastq_path_second  + ' -baseout ' + fastq_prefix + " ILLUMINACLIP:" + adapter_path + ':2:30:10' + ' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 > ' + logfile_fold + '/' + fastq_type + '_trimmomatic.log' + ' 2>&1'
	print cmd_trimmomatic
	os.system(cmd_trimmomatic)


def kallisto_expression(raw_fastq_path_first,raw_fastq_path_second,kallisto_path,kallisto_out_fold,prefix,kallisto_cdna_path,logfile_fold):
	cdna_path_dir = os.path.dirname(kallisto_cdna_path)
	cnd_file_prefix = os.path.splitext(os.path.basename(kallisto_cdna_path))[0]
	kallisto_index_path = cdna_path_dir + '/' + cnd_file_prefix + '.idx'
	cmd_kallisto_index = kallisto_path + " index -i " + kallisto_index_path + ' ' + kallisto_cdna_path + ' > ' +  logfile_fold + '/' + prefix + '_kallisto_index.log' + ' 2>&1'
	cmd_kallisto_quant = kallisto_path + " quant -i " + kallisto_index_path + " -o " + kallisto_out_fold + " " + raw_fastq_path_first + " " + raw_fastq_path_second + ' > ' +  logfile_fold + '/' + prefix + '_kallisto.log' + ' 2>&1'
	print cmd_kallisto_index
	#os.system(cmd_kallisto_index)
	print cmd_kallisto_quant
	os.system(cmd_kallisto_quant)

def hlatyping(raw_fastq_path_first,raw_fastq_path_second,opitype_fold,opitype_out_fold,opitype_ext,prefix):
	os.system("rm -rf %s/*"%opitype_out_fold)
	cmd_hla = 'python ' + opitype_fold + ' -i ' + raw_fastq_path_first + ' ' + raw_fastq_path_second + ' --rna -o ' + opitype_out_fold
	print cmd_hla
	os.system(cmd_hla)
	result_dir=os.listdir(opitype_out_fold)
	print result_dir[0]
	hla_result_path=opitype_out_fold+'/'+result_dir[0]+'/'+result_dir[0]+'_result.tsv'
	print hla_result_path
	cmd_hla_ext = 'python ' + opitype_ext + ' -i ' + hla_result_path + ' -o ' + opitype_out_fold + ' -s ' + prefix
	print cmd_hla_ext
	os.system(cmd_hla_ext)
	print 'hla type process done.'

def mapping_qc_gatk_preprocess(fastq_1_path,fastq_2_path,CPU,STAR_index,alignment_out_fold,prefix,REFERENCE,STAR_path,java_picard_path,GATK_path,dbsnp138,OneKG,mills,GTF_path):
	cmd_STAR=STAR_path + " --runThreadN " + str(CPU) + " --twopassMode Basic --readFilesCommand zcat --readFilesIn " + fastq_1_path + " " +  fastq_2_path + " --genomeDir " + STAR_index + " --outFileNamePrefix " + alignment_out_fold+'/'+prefix + " --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 64424509440 "
	cmd_add_readgroup="java -Xmx4G -jar " + java_picard_path + ' AddOrReplaceReadGroups I=' + alignment_out_fold+'/'+ prefix + 'Aligned.sortedByCoord.out.bam' + ' O=' + alignment_out_fold+'/'+ prefix + '_add.bam' + ' SO=coordinate VALIDATION_STRINGENCY=SILENT RGID=id RGLB=solexa-123 RGPL=illumina RGPU=AXL2342  RGSM=WGC015802 RGCN=bi RGDT=2014-01-20'
	cmd_picard="java -Xmx4G -jar " + java_picard_path + ' MarkDuplicates INPUT=' + alignment_out_fold+'/'+ prefix + '_add.bam' + ' OUTPUT=' + alignment_out_fold+'/'+ prefix + '_add_markdup.bam' + ' METRICS_FILE=' + alignment_out_fold+'/'+prefix + '_dup_qc.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000'
	cmd_buildbamindex="java -Xmx4G -jar " + java_picard_path + ' BuildBamIndex I=' + alignment_out_fold+'/'+ prefix + '_add_markdup.bam' + ' O=' + alignment_out_fold+'/'+ prefix + '_add_markdup.bam.bai' + ' VALIDATION_STRINGENCY=SILENT'
	cmd_SplitNCigarReads="java -Xmx4G -jar " + GATK_path + ' -T SplitNCigarReads -R ' + REFERENCE + ' -I '+ alignment_out_fold+'/'+ prefix + '_add_markdup.bam' + ' -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -o ' + alignment_out_fold+'/'+ prefix + '_SplitNCigarReads.bam'
	cmd_RealignerTargetCreator="java -Xmx4G -jar " + GATK_path + ' -T RealignerTargetCreator -nt 16 -dt NONE -R ' + REFERENCE + ' -I '+ alignment_out_fold+'/'+ prefix + '_SplitNCigarReads.bam' + ' -known ' + OneKG + ' -known ' + mills + ' -filterRNC -o ' + alignment_out_fold+'/'+ prefix + '.intervals'
	cmd_IndelRealigner="java -Xmx4G -jar " + GATK_path + ' -T IndelRealigner -R ' + REFERENCE + ' -I ' + alignment_out_fold+'/'+ prefix + '_SplitNCigarReads.bam' + ' -known ' + OneKG + ' -known ' + mills + ' -targetIntervals ' + alignment_out_fold+'/'+ prefix + '.intervals' + ' -filterRNC -o ' + alignment_out_fold+'/'+ prefix + '_realign.bam'
	cmd_BaseRecalibrator="java -Xmx4G -jar " + GATK_path + ' -T BaseRecalibrator -nct 16 -R ' + REFERENCE + ' -I ' + alignment_out_fold+'/'+ prefix +'_realign.bam' + ' -knownSites ' + OneKG + ' -knownSites ' + mills + ' -knownSites ' + dbsnp138 + ' -o ' + alignment_out_fold + '/' + prefix + '.table'
	cmd_PrintReads="java -Xmx4G -jar " + GATK_path + ' -T PrintReads -nct 16 -dt NONE -R ' + REFERENCE + ' -I ' + alignment_out_fold+'/'+ prefix + '_realign.bam' + ' -BQSR ' + alignment_out_fold + '/' + prefix + '.table' + ' -o ' + alignment_out_fold + '/' + prefix + '_recal.bam'
	print cmd_STAR
	os.system(cmd_STAR)
	print cmd_add_readgroup
	os.system(cmd_add_readgroup)
	print cmd_picard
	os.system(cmd_picard)
	print cmd_buildbamindex
	os.system(cmd_buildbamindex)
	print cmd_SplitNCigarReads
	os.system(cmd_SplitNCigarReads)
	print cmd_RealignerTargetCreator
	os.system(cmd_RealignerTargetCreator)
	print cmd_IndelRealigner
	os.system(cmd_IndelRealigner)
	print cmd_BaseRecalibrator
	os.system(cmd_BaseRecalibrator)
	print cmd_PrintReads
	os.system(cmd_PrintReads)

def GATK_hp(GATK_path,REFERENCE,alignment_out_fold,prefix,CPU,dbsnp138,somatic_out_fold,vcftools_path,vep_path,vep_cache_path,netmhc_out_path,pos_1000G_file):
	cmd_HaplotypeCaller = "java -Xmx4G -jar " + GATK_path + ' -T HaplotypeCaller -nct 32'  + ' -R ' + REFERENCE  + ' -I ' + alignment_out_fold + '/' + prefix + '_'+ 'recal.bam ' + '--dbsnp ' + dbsnp138 + ' -o ' + somatic_out_fold + '/' + prefix + '_'+ 'haplotyper.vcf' + " -dontUseSoftClippedBases -stand_call_conf 20.0"
	print cmd_HaplotypeCaller
	#os.system(cmd_HaplotypeCaller)
	cmd_variantsfilter = "java -Xmx4G -jar " + GATK_path + ' -T VariantFiltration' + ' -R ' + REFERENCE + ' -V ' + somatic_out_fold + '/' + prefix + '_'+ 'haplotyper.vcf' + ' -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ' + somatic_out_fold + '/' + prefix + '_'+ 'haplotyper_filter.vcf'
	#print cmd_variantsfilter
	#os.system(cmd_variantsfilter)	
	cmd_vcftools_snv=vcftools_path + " --vcf " + somatic_out_fold + '/' + prefix + '_'+ 'haplotyper_filter.vcf' + " --remove-filtered-all --remove-indels --recode --recode-INFO-all --minQ 1000 --not-chr chrM --out " + somatic_out_fold + '/' + prefix + '_'+ 'SNVs_only'
	print cmd_vcftools_snv
	os.system(cmd_vcftools_snv)
	cmd_1000G_filter="python ${iTuNES_BIN_PATH}/snv_filter_1000genome.py -i " + somatic_out_fold + '/' + prefix + '_'+ 'SNVs_only.recode.vcf' + ' -g ' + pos_1000G_file + " -o " + somatic_out_fold + ' -s ' + prefix 
	print cmd_1000G_filter
	os.system(cmd_1000G_filter)
	cmd_vcftools_indel=vcftools_path + " --vcf " + somatic_out_fold + '/' + prefix + '_'+ 'haplotyper_filter.vcf' + " --remove-filtered-all --keep-only-indels --recode --recode-INFO-all --not-chr chrM --out " + somatic_out_fold + '/' + prefix + '_'+ 'INDELs_only'
	#print cmd_vcftools_indel
	#os.system(cmd_vcftools_indel)
	cmd_vep=vep_path + " -i " + somatic_out_fold + '/' + prefix + '_'+ 'SNVs_filter_1000.vcf' + " --cache --dir " + vep_cache_path + " --dir_cache " + vep_cache_path + " --force_overwrite --canonical --symbol -o STDOUT --offline | filter_vep --ontology --filter \"CANONICAL is YES and Consequence is missense_variant\" -o " + somatic_out_fold + '/' + prefix + '_'+ 'snv_vep_ann.txt' + " --force_overwrite"
	print cmd_vep
	os.system(cmd_vep)
	cmd_snv="python ${iTuNES_BIN_PATH}/snv2fasta.py -i " + somatic_out_fold + '/' + prefix + '_'+ 'snv_vep_ann.txt' + " -o " + netmhc_out_path + " -s " + prefix
	print cmd_snv
	os.system(cmd_snv)	



def netMHCpan(fasta_file,hla_str,netmhc_out_file,out_dir,split_num,netMHCpan_path,tmp_dir):
	str_proc=r'''
set -x
input_fasta=%s
hla_str=%s
netmhc_out=%s
out_dir=%s
split_num=%s
netMHCpan=%s
tmp=%s
if [ -d ${out_dir}/${tmp} ];then
	rm -rf ${out_dir}/${tmp}
	mkdir -p ${out_dir}/${tmp}
else
	mkdir -p ${out_dir}/${tmp}
fi
if [ -f ${netmhc_out} ];then
	rm ${netmhc_out}
fi
split -l ${split_num} ${input_fasta} ${out_dir}/${tmp}/
filelist=`ls ${out_dir}/${tmp}/`
arr1=(${filelist})
echo ${arr1[@]}
OLD_IFS="$IFS" 
IFS=","
arr2=(${hla_str})
IFS="$OLD_IFS" 
for s in ${arr2[@]}
do
{
	echo $s
	for file_l in ${arr1[@]}
	do
	{
		echo ${file_l}
		$netMHCpan -a $s -f ${out_dir}/${tmp}/${file_l} -l 9 > ${out_dir}/${tmp}/${s}_${file_l}_tmp_netmhc.txt
	} &
	done
	wait
}
done
for file_l in ${arr1[@]}
do
{
	rm ${out_dir}/${tmp}/${file_l}
}
done
filelist1=`ls ${out_dir}/${tmp}/`
for file_r in $filelist1
do
{
	cat ${out_dir}/${tmp}/${file_r} >> ${netmhc_out}
	rm ${out_dir}/${tmp}/${file_r}	
}
done
rm -rf 	${out_dir}/${tmp}
set +x
'''%(fasta_file,hla_str,netmhc_out_file,out_dir,split_num,netMHCpan_path,tmp_dir)
	subprocess.call(str_proc, shell=True, executable='/bin/bash')


def varscan_snv_calling(somatic_mutation_fold,alignment_out_fold,PREFIX,REFERENCE,vep_cache,samtools_path,varscan_path,vep_path,netmhc_out_fold):
	str_proc = r'''
set -e
somatic_mutation=%s
alignment_fold=%s
PREFIX=%s
REFERENCE=%s
vep_cache=%s
netmhc_out=%s
samtools=%s
varscan=%s
vep=%s
if [ ! -d ${somatic_mutation} ];then
	mkdir ${somatic_mutation}	
fi
if [ ! -d ${netmhc_out} ];then
	mkdir ${netmhc_out}	
fi
rm -rf ${somatic_mutation}/*
cd ${somatic_mutation}
mkfifo ${PREFIX}_tumor.fifo
$samtools mpileup -f ${REFERENCE} -q 5 -Q 20 -L 10000 -d 10000 ${alignment_fold}/${PREFIX}_recal.bam > ${PREFIX}_tumor.fifo &
java -jar $varscan mpileup2snp ${PREFIX}_tumor.fifo --output-vcf 1 > ${somatic_mutation}/${PREFIX}_varscan_snv.vcf #--output-vcf 1
rm ${PREFIX}_tumor.fifo
cd ..
$vep -i ${somatic_mutation}/${PREFIX}_varscan_snv.vcf --cache --dir $vep_cache --dir_cache $vep_cache --force_overwrite --canonical --symbol -o STDOUT --offline | filter_vep --ontology --filter "CANONICAL is YES and Consequence is missense_variant" -o ${somatic_mutation}/${PREFIX}_varscan_snv_vep_ann.txt --force_overwrite
python ${iTuNES_BIN_PATH}/snv2fasta.py -i ${somatic_mutation}/${PREFIX}_varscan_snv_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}_varscan
'''%(somatic_mutation_fold,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,samtools_path,varscan_path,vep_path)
	print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')
def varscan_indel_calling(varscan_indel_fold,alignment_out_fold,PREFIX,REFERENCE,vep_cache,samtools_path,varscan_path,vep_path,netmhc_out_fold):
	str_proc = r'''
set -e
varscan_indel=%s
alignment_fold=%s
PREFIX=%s
REFERENCE=%s
vep_cache=%s
netmhc_out=%s
samtools=%s
varscan=%s
vep=%s
if [ ! -d ${varscan_indel} ];then
	mkdir ${varscan_indel}	
fi
if [ ! -d ${netmhc_out} ];then
	mkdir ${netmhc_out}	
fi
#rm -rf ${varscan_indel}/*
#cd ${varscan_indel}
#mkfifo ${PREFIX}_tumor.fifo
#$samtools mpileup -f ${REFERENCE} -q 5 -Q 20 -L 10000 -d 10000 ${alignment_fold}/${PREFIX}_recal.bam  > ${PREFIX}_tumor.fifo &
#java -jar $varscan mpileup2indel ${PREFIX}_tumor.fifo --output-vcf 1 > ${varscan_indel}/${PREFIX}_varscan_indel.vcf #--output-vcf 1
#rm ${PREFIX}_tumor.fifo
#cd ..
$vep -i ${varscan_indel}/${PREFIX}_varscan_indel.vcf --cache --dir $vep_cache --dir_cache $vep_cache --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is coding_sequence_variant" -o ${varscan_indel}/${PREFIX}_varscan_indel_vep_ann.txt --force_overwrite
python ${iTuNES_BIN_PATH}/varscandel2fasta.py -i ${varscan_indel}/${PREFIX}_varscan_indel_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}_varscan
python ${iTuNES_BIN_PATH}/varscanins2fasta.py -i ${varscan_indel}/${PREFIX}_varscan_indel_vep_ann.txt  -o ${netmhc_out} -s ${PREFIX}_varscan
'''%(varscan_indel_fold,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,samtools_path,varscan_path,vep_path)
	print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')

def strelka_indel_calling(strelka_out_fold,strelka_path,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,CPU,vep_path):
	str_proc2=r'''
set -e
strelka_fold=%s
strelka_path=%s
alignment_fold=%s
PREFIX=%s
REFERENCE=%s
vep_cache=%s
netmhc_out=%s
cpu=%s
vep=%s
if [ -d ${strelka_fold} ];then
	rm -rf ${strelka_fold}
fi
python ${strelka_path}/configureStrelkaGermlineWorkflow.py --bam=${alignment_fold}/${PREFIX}_recal.bam --referenceFasta=${REFERENCE} --config=${strelka_path}/configureStrelkaGermlineWorkflow.py.ini --runDir=${strelka_fold} --rna
python ${strelka_fold}/runWorkflow.py -m local -j $cpu -q ${PREFIX}_strelka -g 32 --quiet
$vep -i ${strelka_fold}/results/variants/germline.indels.vcf.gz --cache --dir ${vep_cache} --dir_cache ${vep_cache} --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is coding_sequence_variant" -o ${strelka_fold}/${PREFIX}_strelka_indel_vep_ann.txt --force_overwrite
python ${iTuNES_BIN_PATH}/varscandel2fasta.py -i ${strelka_fold}/${PREFIX}_strelka_indel_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}_strelka
python ${iTuNES_BIN_PATH}/varscanins2fasta.py -i ${strelka_fold}/${PREFIX}_strelka_indel_vep_ann.txt  -o ${netmhc_out} -s ${PREFIX}_strelka
'''%(strelka_out_fold,strelka_path,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,CPU,vep_path)
	print str_proc2
	subprocess.call(str_proc2, shell=True, executable='/bin/bash')
	
def varscan_neo(snv_fasta_file,hla_str,snv_netmhc_out_file,netmhc_out_fold,split_num,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_fold,netMHCpan_path):
	netMHCpan(snv_fasta_file,hla_str,snv_netmhc_out_file,netmhc_out_fold,split_num,netMHCpan_path,'tmp_snv')
	str_proc1=r'''
PREFIX=%s
netmhc_out=%s
Exp_file=%s
Binding_Aff_Fc_Cutoff=%d
Binding_Aff_Cutoff=%d
Fpkm_Cutoff=%d
hla_str=%s
netctl_fold=%s
python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i ${netmhc_out}/${PREFIX}_snv_netmhc.txt -g ${netmhc_out}/${PREFIX}_snv.fasta -o ${netmhc_out} -s ${PREFIX}_snv -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
python ${iTuNES_BIN_PATH}/netCTLPAN.py -i ${netmhc_out}/${PREFIX}_snv_final_neo_candidate.txt -o ${netctl_fold} -s ${PREFIX}_snv
'''%(prefix,netmhc_out_fold,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,hla_str,netctl_fold)
	print str_proc1
	subprocess.call(str_proc1, shell=True, executable='/bin/bash')

def indel_neo(PREFIX,netmhc_out_fold,indel_fasta_file,hla_str,indel_netmhc_out_file,split_num,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_fold,netMHCpan_path):
	str_proc3=r'''
PREFIX=%s
netmhc_out=%s
#cat ${netmhc_out}/${PREFIX}_strelka_del.fasta > ${netmhc_out}/${PREFIX}_indel.fasta
#cat ${netmhc_out}/${PREFIX}_strelka_ins.fasta >> ${netmhc_out}/${PREFIX}_indel.fasta
cat ${netmhc_out}/${PREFIX}_gatk_del.fasta >> ${netmhc_out}/${PREFIX}_indel.fasta
cat ${netmhc_out}/${PREFIX}_gatk_ins.fasta >> ${netmhc_out}/${PREFIX}_indel.fasta
'''%(PREFIX,netmhc_out_fold)
	subprocess.call(str_proc3, shell=True, executable='/bin/bash')
	netMHCpan(indel_fasta_file,hla_str,indel_netmhc_out_file,netmhc_out_fold,split_num,netMHCpan_path,'tmp_indel')
	str_proc4=r'''
set -e
PREFIX=%s
netmhc_out=%s
Exp_file=%s
Binding_Aff_Fc_Cutoff=%s
Binding_Aff_Cutoff=%s
Fpkm_Cutoff=%s
hla_str=%s
netctl_fold=%s
python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i ${netmhc_out}/${PREFIX}_indel_netmhc.txt -g ${netmhc_out}/${PREFIX}_indel.fasta -o ${netmhc_out} -s ${PREFIX}_indel -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
python ${iTuNES_BIN_PATH}/netCTLPAN.py -i ${netmhc_out}/${PREFIX}_indel_final_neo_candidate.txt -o ${netctl_fold} -s ${PREFIX}_indel
'''%(PREFIX,netmhc_out_fold,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,hla_str,netctl_fold)
	#print str_proc4
	subprocess.call(str_proc4, shell=True, executable='/bin/bash')

def gene_fusion(fastq_1_path,fastq_2_path,prefix,gene_fusion_fold_path,eric_db_path,eric_path,cpu,netmhc_out_fold,genefuison_fasta_file,hla_str,genefuison_netmhc_out_file,split_num,netMHCpan_path,Binding_Aff_Cutoff,Fpkm_Cutoff):
	str_proc='''
T_fastq_1=%s
T_fastq_2=%s
prefix=%s
gene_fusion_fold=%s
eric_db=%s
eric_path=%s
cpu=%d
netmhc_out=%s
$eric_path -db $eric_db --refid homo_sapiens -name ${prefix}_genefusion -o result $T_fastq_1 $T_fastq_2 -p $cpu
cd ..
python ${iTuNES_BIN_PATH}/genefusion2fasta.py -i ${gene_fusion_fold}/result/${prefix}_genefusion.results.total.tsv -c ${GeneFusion_Cutoff} -o ${netmhc_out} -s ${prefix}
'''%(fastq_1_path,fastq_2_path,prefix,gene_fusion_fold_path,eric_db_path,eric_path,cpu,netmhc_out_fold)
	netMHCpan(genefuison_fasta_file,hla_str,genefuison_netmhc_out_file,netmhc_out_fold,split_num,netMHCpan_path,'tmp_genefuison')
	str_proc1='''
netmhc_out=%s
prefix=%s
Binding_Aff_Cutoff=%s
Fpkm_Cutoff=%s
hla_str=%s
python ${iTuNES_BIN_PATH}/gf_netMHC_result_parse.py -i ${netmhc_out}/${prefix}_genefusion_netmhc.txt -t ${netmhc_out}/${prefix}_gene_fusion.fasta -g ${gene_fusion_fold}/result/${prefix}_genefusion.results.total.tsv -o ${netmhc_out} -s ${prefix} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
python ${iTuNES_BIN_PATH}/gf_netctl.py -i netmhc/${prefix}_genefusion_final_neo_candidate.txt -o netctl -s ${prefix}_gf
'''%(netmhc_out_fold,prefix,Binding_Aff_Cutoff,Fpkm_Cutoff,hla_str)

def plot_results(X, Y_, means, covariances, index, title, gmm_classification_file):
	color_iter = itertools.cycle(['navy', 'c'])
	splot = plt.subplot(1, 1, 1 + index)
	#splot=plt.plot()
	for i, (mean, covar, color) in enumerate(zip(means, covariances, color_iter)):
		print i,color
		v, w = linalg.eigh(covar)
		v = 2. * np.sqrt(2.) * np.sqrt(v)
		u = w[0] / linalg.norm(w[0])
		# as the DP will not use every component it has access to
		# unless it needs it, we shouldn't plot the redundant
		# components.
		if not np.any(Y_ == i):
			continue
		plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color)
		# Plot an ellipse to show the Gaussian component
		angle = np.arctan(u[1] / u[0])
		angle = 180. * angle / np.pi  # convert to degrees
		ell = mpl.patches.Ellipse(mean, v[0], v[1], 180. + angle, color=color)
		ell.set_clip_box(splot.bbox)
		ell.set_alpha(0.5)
		splot.add_artist(ell)
	plt.xlim(-1.0, 1.8)
	plt.ylim(-.5, .8)
	plt.xlabel('first component')
	plt.ylabel('second component')
	#plt.xticks()
	#plt.yticks()
	plt.title(title)
	plt.savefig(gmm_classification_file)

def immunogenicity_score_calculate(final_neo_file,gmm_classification_file,immunogenicity_score_ranking,immunogenicity_gmm_score_ranking):
	data_neo = pd.read_table(final_neo_file,header=0,sep='\t')
	neoantigen_infor = []
	for i in range(len(data_neo.Gene)):
		neo_infor = data_neo["Gene"][i]+'_'+data_neo["AA_change"][i]+'_'+data_neo["MT_pep"][i]+'_'+data_neo["WT_pep"][i]
		neoantigen_infor.append(neo_infor)
	data_neo["neoantigen_infor"] = neoantigen_infor

	f_affinity_rank_wt=lambda x:1-(1/(1+math.pow(math.e,5*(x-2))))/2
	f_affinity_rank_mt=lambda x:1/(1+math.pow(math.e,5*(x-2)))
	k=1
	f_TPM=lambda x:math.tanh(x/k)
	f_normal_TPM=lambda x:1-math.tanh(x/k)
	aff_mt_rank_score=data_neo.MT_Binding_level.apply(f_affinity_rank_mt)
	aff_wt_rank_score=data_neo.WT_Binding_level.apply(f_affinity_rank_wt)
	#allele_frequency_score=data_neo.variant_allele_frequency
	netchop_score=data_neo.combined_prediction_score
	#cellular_prevalence_score=data_neo.cellular_prevalence
	tpm_score=data_neo.tpm.apply(f_TPM)
	#tpm_normal_score=data_neo.iloc[:,22].apply(f_normal_TPM)
	immunogenicity_score=[]
	for i in range(len(aff_mt_rank_score)):
		IS=aff_mt_rank_score[i]*aff_wt_rank_score[i]*tpm_score[i]*netchop_score[i]
		immunogenicity_score.append(IS)
	data_feature_select=pd.DataFrame()
	data_feature_select["neoantigen_infor"]=neoantigen_infor
	data_feature_select["aff_mt_rank_score"]=aff_mt_rank_score
	data_feature_select["aff_wt_rank_score"]=aff_wt_rank_score
	data_feature_select["tpm_score"]=tpm_score
	#data_feature_select["tpm_normal_score"]=tpm_normal_score
	#data_feature_select["allele_frequency_score"]=allele_frequency_score
	data_feature_select["netchop_score"]=netchop_score
	#data_feature_select["cellular_prevalence_score"]=cellular_prevalence_score
	data_feature_select["immunogenicity_score"]=immunogenicity_score
	data_feature_select_sorted=data_feature_select.sort_values(["immunogenicity_score"],ascending=False)
	data_feature_select_sorted.to_csv(immunogenicity_score_ranking,header=1,index=0,sep='\t')
	X_neo = data_feature_select.values[:,1:5]
	####pca on non-standardized data of all 4 feature
	pca = PCA(n_components=2).fit(X_neo)
	X_neo_pca = pca.transform(X_neo)
	#print "PCA result on non-standardized data"
	#print "explained_variance_ratio",pca.explained_variance_ratio_
	#print "explained_variance",pca.explained_variance_
	#print pca.n_components_
	#print pca.components_
	###plot GMM ellipse
	# Fit a Gaussian mixture with EM using five components
	gmm = mixture.GaussianMixture(n_components=2, covariance_type='full',n_init=5).fit(X_neo_pca)
	plot_results(X_neo_pca, gmm.predict(X_neo_pca), gmm.means_, gmm.covariances_, 0,'Gaussian Mixture',gmm_classification_file)
	#plot_results(X_neo_pca, gmm.predict(X_neo_pca), gmm.means_, gmm.covars_, 0, 'Gaussian Mixture')
	predict_label=gmm.predict(X_neo_pca)
	predict_prob=gmm.predict_proba(X_neo_pca)
	predict_positive_prob=[]
	for i in range(len(predict_prob)):
		pos_prob=predict_prob[i][1]
		predict_positive_prob.append(pos_prob)
	data_feature_select["gmm_label"]=predict_label
	data_gmm_filter_1=data_feature_select[data_feature_select["gmm_label"]==1]
	data_gmm_filter_1_scoreAve=data_gmm_filter_1.immunogenicity_score.sum()/len(data_gmm_filter_1.immunogenicity_score)
	data_gmm_filter_0=data_feature_select[data_feature_select["gmm_label"]==0]
	data_gmm_filter_0_scoreAve=data_gmm_filter_0.immunogenicity_score.sum()/len(data_gmm_filter_0.immunogenicity_score)
	if data_gmm_filter_0_scoreAve > data_gmm_filter_1_scoreAve:
		data_gmm_filter_0.gmm_label=1
		data_gmm_filter_1.gmm_label=0
	else:
		pass
	data_labeled=pd.concat([data_gmm_filter_1,data_gmm_filter_0])
	data_label_sorted=data_labeled.sort_values(["immunogenicity_score"],ascending=False)
	#data_label_sorted.to_csv("gmm_score_sorted.txt",header=1,index=0,sep='\t')
	data_gmm_positive=data_label_sorted[data_label_sorted["gmm_label"]==1]
	data_gmm_positive.to_csv(immunogenicity_gmm_score_ranking,header=1,index=0,sep='\t')


