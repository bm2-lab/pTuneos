import subprocess
import os,sys,time

def hg_liftover(crossmap_path,liftover_file,vcf_file,reference,vcf_fold,vcftools_path,prefix):
	str_proc=r'''
crossmap=%s
chain_file=%s
vcf_file=%s
reference=%s
vcf_fold=%s
vcftools=%s
PREFIX=%s
$crossmap vcf ${chain_file} ${vcf_file} $reference ${vcf_fold}/${PREFIX}_liftover38.vcf
$vcftools --vcf ${vcf_fold}/${PREFIX}_liftover38.vcf --remove-indels --recode --recode-INFO-all --not-chr chrM --out ${vcf_fold}/${PREFIX}_SNVs_only
$vcftools --vcf ${vcf_fold}/${PREFIX}_liftover38.vcf --keep-only-indels --recode --recode-INFO-all --not-chr chrM --out ${vcf_fold}/${PREFIX}_INDELs_only
'''%(crossmap_path,liftover_file,vcf_file,reference,vcf_fold,vcftools_path,prefix)
	print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')
def hg_no_liftover(vcf_file,vcf_fold,vcftools_path,prefix):
	str_proc=r'''
vcf_file=%s
vcf_fold=%s
vcftools=%s
PREFIX=%s
$vcftools --vcf ${vcf_file} --remove-filtered-all --remove-indels --recode --recode-INFO-all --not-chr chrM --out ${vcf_fold}/${PREFIX}_SNVs_only
$vcftools --vcf ${vcf_file} --remove-filtered-all --keep-only-indels --recode --recode-INFO-all --not-chr chrM --out ${vcf_fold}/${PREFIX}_INDELs_only
'''%(vcf_file,vcf_fold,vcftools_path,prefix)
	print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')



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
		$netMHCpan -a $s -f ${out_dir}/${tmp}/${file_l} -l 8,9,10,11 > ${out_dir}/${tmp}/${s}_${file_l}_tmp_netmhc.txt
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



def snv_neo(vep_path,snv_vcf_file,vep_cache,vcf_fold,netmhc_out_fold,snv_fasta_file,netmhc_out_file,split_num,netMHCpan_path,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,hla_str,netctl_fold):
	str_proc1=r'''
PREFIX=%s
vep=%s
snv_vcf_file=%s
VEP_CACHE=%s
vcf_fold=%s
netmhc_out=%s
$vep -i ${snv_vcf_file} --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep.pl --ontology --filter "Consequence is missense_variant" -o ${vcf_fold}/${PREFIX}_snv_vep_ann.txt --force_overwrite
python ${iTuNES_BIN_PATH}/snv2fasta.py -i ${vcf_fold}/${PREFIX}_snv_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}
'''%(prefix,vep_path,snv_vcf_file,vep_cache,vcf_fold,netmhc_out_fold)
	print str_proc1
	#subprocess.call(str_proc1, shell=True, executable='/bin/bash')
	print "start netmhcpan"
	#netMHCpan(snv_fasta_file,hla_str,netmhc_out_file,netmhc_out_fold,split_num,netMHCpan_path,"tmp_snv")
	str_proc2=r'''
PREFIX=%s
netmhc_out=%s
Exp_file=%s
Binding_Aff_Fc_Cutoff=%d
Binding_Aff_Cutoff=%d
Fpkm_Cutoff=%d
hla_str=%s
netctl_fold=%s
#python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i ${netmhc_out}/${PREFIX}_snv_netmhc.txt -g ${netmhc_out}/${PREFIX}_snv.fasta -o ${netmhc_out} -s ${PREFIX}_snv -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
python ${iTuNES_BIN_PATH}/netCTLPAN.py -i ${netmhc_out}/${PREFIX}_snv_final_neo_candidate.txt -o ${netctl_fold} -s ${PREFIX}
#python ${iTuNES_BIN_PATH}/neoantigens_annotation.py -i ${netctl_fold}/${PREFIX}_snv_netctl_concact.txt -o annotation -s ${PREFIX}
'''%(prefix,netmhc_out_fold,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,hla_str,netctl_fold)
	subprocess.call(str_proc2, shell=True, executable='/bin/bash')


def indel_neo(vep_path,indel_vcf_file,vep_cache,vcf_fold,netmhc_out_fold,indel_fasta_file,hla_str,netmhc_out_file,split_num,netMHCpan_path,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_fold):
	str_proc1=r'''
PREFIX=%s
vep=%s
indel_vcf_file=%s
VEP_CACHE=%s
vcf_fold=%s
netmhc_out=%s
$vep -i ${indel_vcf_file} --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep.pl --ontology --filter "Consequence is coding_sequence_variant" -o ${vcf_fold}/${PREFIX}_indel_vep_ann.txt --force_overwrite
python ${iTuNES_BIN_PATH}/varscandel2fasta.py -i ${vcf_fold}/${PREFIX}_indel_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}
python ${iTuNES_BIN_PATH}/varscanins2fasta.py -i ${vcf_fold}/${PREFIX}_indel_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}
cat ${netmhc_out}/${PREFIX}_del.fasta > ${netmhc_out}/${PREFIX}_indel.fasta
cat ${netmhc_out}/${PREFIX}_ins.fasta >> ${netmhc_out}/${PREFIX}_indel.fasta
'''%(prefix,vep_path,indel_vcf_file,vep_cache,vcf_fold,netmhc_out_fold)
	#subprocess.call(str_proc1, shell=True, executable='/bin/bash')
	netMHCpan(indel_fasta_file,hla_str,netmhc_out_file,netmhc_out_fold,split_num,netMHCpan_path,"tmp_indel")
	str_proc2=r'''
PREFIX=%s
netmhc_out=%s
Exp_file=%s
Binding_Aff_Fc_Cutoff=%d
Binding_Aff_Cutoff=%d
Fpkm_Cutoff=%d
hla_str=%s
netctl_fold=%s
python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i ${netmhc_out}/${PREFIX}_indel_netmhc.txt -g ${netmhc_out}/${PREFIX}_indel.fasta -o ${netmhc_out} -s ${PREFIX} -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff}
python ${iTuNES_BIN_PATH}/netCTLPAN.py -i ${netmhc_out}/${PREFIX}_indel_final_neo_candidate.txt -o ${netctl_fold} -s ${PREFIX}
'''%(prefix,netmhc_out_fold,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,hla_str,netctl_fold)	
	#subprocess.call(str_proc2, shell=True, executable='/bin/bash')















