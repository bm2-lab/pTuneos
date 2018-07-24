import argparse

def ParseVCF(p_vcf):
	p_vcf_inp = p_vcf.add_argument_group('Input options')
	p_vcf_inp.add_argument('-i', '--config_file', dest='Config_file', required=True, help='input configure file of programme(required)')