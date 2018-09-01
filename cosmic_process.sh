help_info(){
	echo "cosmic_process.sh: A tool to format the COSMIC VCF file compatible for Mutect2"
	echo " 	Usage: bash cosmic_process.sh <-i CosmicCodingMuts.vcf> <-o FormattedCosmicMuts.vcf> <-d Homo_sapiens_assembly38.dic> "
	echo "	-i  Input CosmicCodingMuts.vcf file as downloaded from COSMIC "
	echo "	-o	Output formatted COSMIC VCF file"
	echo "	-p  your path to picard"
	echo "	-d	Path to the GATK hg19 bundle's sequence dictionary file"

}


if [ $# -lt 2 ];then
	help_info
	exit 1
fi

TEMP=`getopt -o i:o:p:d: \
--long input_vcf:,output_vcf:,picard_path:,dic_file: \
	-n 'cosmic_process.sh' -- "$@"`

if [ $? != 0 ];then
	echo "Terminating..." >&2 ; exit 1 ; fi


eval set -- "$TEMP"
while true
do
	case "$1" in
		-i | --input_vcf) input_comisc_vcf=$2;shift 2;;
		-o | --output_vcf) output_cosmic_vcf=$2;shift 2;;
		-p | --picard_path) PICARD=$2;shift 2;;
		-d | --dic) dic_file=$2;shift 2;;
		--) shift;break;;
		*) echo "Internal error!";exit 1;;
	esac

done

awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' $input_comisc_vcf | sed 's/chrMT/chrM/g' - > CosmicCodingMuts_chr_M.vcf

java -jar $PICARD SortVcf I=CosmicCodingMuts_chr_M.vcf O=$output_cosmic_vcf SEQUENCE_DICTIONARY=$dic_file

rm CosmicCodingMuts_chr_M.vcf

rm ${output_cosmic_vcf}.idx

bgzip ${output_cosmic_vcf} > ${output_cosmic_vcf}.gz
tabix ${output_cosmic_vcf}.gz
