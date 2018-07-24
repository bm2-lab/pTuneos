import argparse

def ParsePAIRNOMATCHDNA(p_penmdna):
	p_penmdna_inp = p_penmdna.add_argument_group('Input options')
	p_penmdna_inp.add_argument('-i', '--config_file', dest='config_file', required=True, help='Input configure file for tumor_only neoantigens identification(required)')
