import argparse
import os
import sys
from src import ui
import shutil

def main():
	
	str_prog = 'iTuNES'
	str_usage = 'python iTuNES.py <command> [options]'
	str_desc = r'''
some description
'''

	p = argparse.ArgumentParser(prog=str_prog,
                                 usage=str_usage,
                                 description=str_desc,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,add_help=True)

	p.add_argument('-v', '--version', action='version', version='%(prog)s 0.3.2(dev)')

	sp = p.add_subparsers(title='Command', metavar='')

	p_PairMatchDna = sp.add_parser('PairMatchDna',
						description='detecting tumor neoantigens using matched pair-end WGS/WES sequencing data',
						usage='python iTuNES.py pairmatchdna [options]',
						help='detecting tumor neoantigens using pair-end WGS/WES sequencing data')
	ui.ParsePAIRMATCHDNA(p_PairMatchDna)
	p_PairNoMatchDna = sp.add_parser('PairNomatchDna',
                         description='detecting tumor neoantigens using pair-end WGS/WES sequencing data',
                         usage='python iTuNES.py pairmatchdna [options]',
                         help='detecting tumor neoantigens using pair-end WGS/WES sequencing data')
	ui.ParsePAIRNOMATCHDNA(p_PairNoMatchDna)
	p_GeneFusion = sp.add_parser('GeneFusion',
                         description='detecting tumor neoantigens using single-end WGS/WES sequencing data',
                         usage='python iTuNES.py GeneFusion [options]',
                         help='detecting tumor neoantigens using single-end WGS/WES sequencing data')
	ui.ParseGENEFUSION(p_GeneFusion)
	p_PairRna = sp.add_parser('PairRna',
                         description='detecting tumor neoantigens using pair-end RNA sequencing data',
                         usage='python iTuNES.py pairmatchdna [options]',
                         help='detecting tumor neoantigens using pair-end RNA sequencing data')
	ui.ParsePAIRRNA(p_PairRna)
	p_VCF = sp.add_parser('VCF',
                         description='detecting tumor neoantigens using somatic mutaion data in VCF file format',
                         usage='python iTuNES.py vcf [options]',
                         help='detecting tumor neoantigens using somatic mutaion data in VCF file format')
	ui.ParseVCF(p_VCF)
	if len(sys.argv) == 1:
		p.print_help()
		sys.exit(1)
	elif len(sys.argv) == 2:
		if sys.argv[1] == '-h' or sys.argv[1] == '--help':
           		p.print_help()
		if sys.argv[1] == 'PairMatchDna':
			p_PairMatchDna.print_help()
			sys.exit(1)
		if sys.argv[1] == 'PairNomatchDna':
			p_PairNoMatchDna.print_help()
			sys.exit(1)
		if sys.argv[1] == 'p_GeneFusion':
			p_GeneFusion.print_help()
			sys.exit(1)
		if sys.argv[1] == 'PairRna':
			p_PairRna.print_help()
			sys.exit(1)
		if sys.argv[1] == 'VCF':
			p_VCF.print_help()
			sys.exit(1)
	opts = p.parse_args()
	if sys.argv[1] == 'PairMatchDna':
		from src.core import pairendMDNA 
		pairendMDNA.PEMD(opts)
	if sys.argv[1] == 'PairNomatchDna':
		from src.core import pairendNMDNA 
		pairendNMDNA.PENMD(opts)
	if sys.argv[1] == 'GeneFusion':
		from src.core import GeneFusionNeo 
		GeneFusionNeo.GF(opts)
	if sys.argv[1] == 'PairRna':
		from src.core import pairendRNA 
		pairendRNA.PERNA(opts)
	if sys.argv[1] == 'VCF':
		from src.core import VariantCallingFormat 
		VariantCallingFormat.Vcf(opts)

if __name__ == '__main__':
	main()
