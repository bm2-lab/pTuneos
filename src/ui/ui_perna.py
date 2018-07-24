import argparse

def ParsePAIRRNA(p_perna):
	p_perna_inp = p_perna.add_argument_group('Input options')
	p_perna_inp.add_argument('-i', '--config_file', dest='Config_file', required=True, help='input configure file of programme(required)')

