import subprocess
import os,sys,time

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

def gene_fusion(fastq_1_path,fastq_2_path,prefix,gene_fusion_fold_path,eric_db_path,eric_path,cpu,netmhc_out_fold,genefuison_fasta_file,hla_str,genefuison_netmhc_out_file,split_num,netMHCpan_path,Binding_Aff_Cutoff,Fpkm_Cutoff):
	str_proc='''
T_fastq_1=%s
T_fastq_2=%s
prefix=%s
gene_fusion_fold=%s
eric_db=%s
eric_path=%s
cpu=%s
netmhc_out=%s
$eric_path -db $eric_db --refid homo_sapiens -name ${prefix}_genefusion -o ${gene_fusion_fold}/result $T_fastq_1 $T_fastq_2 -p $cpu
python ${iTuNES_BIN_PATH}/genefusion2fasta.py -i ${gene_fusion_fold}/result/${prefix}_genefusion.results.total.tsv -c ${GeneFusion_Cutoff} -o ${netmhc_out} -s ${prefix}
'''%(fastq_1_path,fastq_2_path,prefix,gene_fusion_fold_path,eric_db_path,eric_path,cpu,netmhc_out_fold)
	print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')
	#netMHCpan(genefuison_fasta_file,hla_str,genefuison_netmhc_out_file,netmhc_out_fold,split_num,netMHCpan_path,'tmp_genefuison')
	str_proc1='''
netmhc_out=%s
prefix=%s
Binding_Aff_Cutoff=%s
Fpkm_Cutoff=%s
hla_str=%s
python ${iTuNES_BIN_PATH}/gf_netMHC_result_parse.py -i ${netmhc_out}/${prefix}_genefusion_netmhc.txt -t ${netmhc_out}/${prefix}_gene_fusion.fasta -g ${gene_fusion_fold}/result/${prefix}_genefusion.results.total.tsv -o ${netmhc_out} -s ${prefix} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
python ${iTuNES_BIN_PATH}/gf_netctl.py -i netmhc/${prefix}_genefusion_final_neo_candidate.txt -o netctl -s ${prefix}_gf
'''%(netmhc_out_fold,prefix,Binding_Aff_Cutoff,Fpkm_Cutoff,hla_str)
	print str_proc1
	#subprocess.call(str_proc1, shell=True, executable='/bin/bash')










