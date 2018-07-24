import argparse

def ParsePAIRMATCHDNA(p_pedna):
	p_pedna_inp = p_pedna.add_argument_group('Input options')
	p_pedna_inp.add_argument('-i', '--config_file', dest='Config_file', required=True, help='input configure file of programme(required)')	
	
