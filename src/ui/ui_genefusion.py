import argparse

def ParseGENEFUSION(p_gf):
	p_gf_inp = p_gf.add_argument_group('Input options')
	p_gf_inp.add_argument('-i', '--config_file', dest='Config_file', required=True, help='input configure file of programme(required)')
