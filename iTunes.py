import argparse
import os
import sys
from src import ui
import shutil

def main():
	
	str_prog = 'iTuNES'
	str_usage = 'python iTuNES.py <command> [options]'
	str_desc = r'''
identification of personalized Tumor neoantigens from next-generation sequencing data
'''

	p = argparse.ArgumentParser(prog=str_prog,
                                 usage=str_usage,
                                 description=str_desc,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,add_help=True)

	p.add_argument('-v', '--version', action='version', version='%(prog)s 1.0.0(dev)')

	sp = p.add_subparsers(title='Command', metavar='')

	p_PairMatchDna = sp.add_parser('WES',
						description='detecting tumor neoantigens using matched pair-end WGS/WES sequencing data',
						usage='python iTuNES.py WES [options]',
						help='detecting tumor neoantigens using pair-end WGS/WES sequencing data')
	ui.ParsePAIRMATCHDNA(p_PairMatchDna)
	p_VCF = sp.add_parser('VCF',
                         description='detecting tumor neoantigens using somatic mutaion data in VCF file format',
                         usage='python iTuNES.py VCF [options]',
                         help='detecting tumor neoantigens using somatic mutaion data in VCF file format')
	ui.ParseVCF(p_VCF)
	if len(sys.argv) == 1:
		p.print_help()
		sys.exit(1)
	elif len(sys.argv) == 2:
		if sys.argv[1] == '-h' or sys.argv[1] == '--help':
           		p.print_help()
		if sys.argv[1] == 'WES':
			p_PairMatchDna.print_help()
			sys.exit(1)
		if sys.argv[1] == 'VCF':
			p_VCF.print_help()
			sys.exit(1)
	opts = p.parse_args()
	if sys.argv[1] == 'WES':
		from src.core import pairendMDNA 
		pairendMDNA.PEMD(opts)
	if sys.argv[1] == 'VCF':
		from src.core import VariantCallingFormat 
		VariantCallingFormat.Vcf(opts)

if __name__ == '__main__':
	main()
