#!/usr/bin/env python

'''
Created on 06.16.2016
@author: Dorjee Tamang
'''

import os
from optparse import OptionParser
import re
import logging
logging.basicConfig(level=logging.INFO)

from netchop_3_1_executable import predict as netchop_predict
from netctl_1_1_executable import predict as netctl_predict
from netctlpan_1_1_executable import predict as netctlpan_predict
from allele_info import NetCTLpanAlleleData

def is_num(s):
    try:
        float(s)
        return True
    except:
        return False
            
class PredictNetChop(object):
    
    def __init__(self):
        self.threshold = None
             
    def predict(self, options, args):
        """Returns the prediction result."""        
        self.args_num_validation(len(args))
        input_fasta_list = self.fasta_to_list(args[0])
        self.input_fasta_validation(input_fasta_list)
        
        scores = {}
        if options.method:
            plot_name = options.plot_name or "plot"
            if options.method == 'netchop':
                self.threshold = options.threshold if options.threshold else '0.500000'
                pred_method = options.pred_method if options.pred_method else '0'
                self.netchop_validation(threshold=self.threshold, method=pred_method)
                scores = netchop_predict(input_fasta_list, threshold=self.threshold, method=pred_method)
            elif options.method == 'netctl':
                self.threshold = options.threshold if options.threshold else '0.750000'
                tap_weight = options.tap_weight if options.tap_weight else '0.050000'
                cleavage_weight = options.cleavage_weight if options.cleavage_weight else '0.100000'
                supertype = options.supertype if options.supertype else 'A1'
                self.netctl_validation(tap_weight=tap_weight, cleavage_weight=cleavage_weight, supertype=supertype, threshold=self.threshold)
                scores = netctl_predict(input_fasta_list, tap_weight=tap_weight, cleavage_weight=cleavage_weight, supertype=supertype, threshold=self.threshold)
            elif options.method == 'netctlpan':
                self.threshold = options.threshold if options.threshold else '-99.900002'
                tap_weight = options.tap_weight if options.tap_weight else '0.025000'
                cleavage_weight = options.cleavage_weight if options.cleavage_weight else '0.225000'
                epitope_threshold = options.epitope_threshold if options.epitope_threshold else '1.000000'
                allele = options.allele if options.allele else 'HLA-A02:01'
                length = options.length if options.length else '9'
                self.netctlpan_validation(allele, length, tap_weight, cleavage_weight, epitope_threshold, self.threshold)
                scores = netctlpan_predict(input_fasta_list, allele=allele, length=length, tap_weight=tap_weight, cleavage_weight=cleavage_weight, epitope_threshold=epitope_threshold, threshold=self.threshold)
        else:
            raise ValueError("A method name must be provided. Please check your input.")
        
        logging.debug(scores)
        self.print_results(scores, options.method)
        
        # generate a plot if 'noplot' option is not invoked
        if options.noplot == False:
            self.create_png(scores, options.method, self.threshold, plot_name=plot_name)
            print '* A plot: "%s.png" has been generated.' % plot_name
    
    def input_fasta_validation(self, input_fasta_list):
        if len(input_fasta_list)<1:
            raise ValueError("A protein sequence must be provided. Please check your data file.")
        if len(input_fasta_list)>1:
            raise ValueError("NetChop Standalone accept one protein sequence per time. Please check your input.")
        
    def netchop_validation(self, threshold, method):
        if not method in ['0', '1', '20S']:
            raise ValueError("pred_method must be method id: '0', '1', corresponding method name: '20S', 'C term'.")
        if not is_num(threshold):        
            raise ValueError("input threshold must be number.")            

    
    def netctl_validation(self, tap_weight, cleavage_weight, supertype, threshold):
        if not is_num(tap_weight):        
            raise ValueError("input tap_weight must be number.")
        if not is_num(cleavage_weight):        
            raise ValueError("input cleavage_weight must be number.")
        if not is_num(threshold):        
            raise ValueError("input threshold must be number.")
        if not supertype in ['A1', 'A2', 'A3', 'A24', 'A26', 'B7', 'B8', 'B27', 'B39', 'B44', 'B58', 'B62']:
            raise ValueError("input supertype must be one of these: %s." % ', '.join(['A1', 'A2', 'A3', 'A24', 'A26', 'B7', 'B8', 'B27', 'B39', 'B44', 'B58', 'B62']))
        
    def netctlpan_validation(self, allele, length, tap_weight, cleavage_weight, epitope_threshold, threshold):
    
        if not is_num(tap_weight):        
            raise ValueError("input tap_weight must be number.")
        if not is_num(cleavage_weight):        
            raise ValueError("input cleavage_weight must be number.")
        if not is_num(epitope_threshold):        
            raise ValueError("input epitope_threshold must be number.")
        if not is_num(threshold):        
            raise ValueError("input threshold must be number.")
        ctlpanad = NetCTLpanAlleleData()        
        allele_list = [alle for species in ctlpanad.get_species_list() for alle in ctlpanad.get_allele_names_for_species(species)]        
        if not allele in allele_list:
            raise ValueError("input allele must be one of these allele:\n%s." % '\n'.join(allele_list))
        length_list = ctlpanad.get_allowed_peptide_lengths(allele)
        if not int(length) in length_list:
            raise ValueError("input length must be one of these: %s." % (length_list))
    
    def fasta_to_list(self, fasta_file):
        fasta_list = []
        with open(fasta_file, "r") as infile:
            file_content = infile.read()
            if '>' not in file_content:
                raise ValueError('Expected ">" not found. Input file must be in fasta format.')

            sequences = file_content.split('>')
            for s_raw in sequences:
                s = s_raw.strip()
                if len(s) == 0: continue
                end_of_name = s.find('\n')                
                sequence_name = s[:end_of_name].strip()
                if not sequence_name:
                    raise ValueError("No sequence_name. Please check your fasta file.")
                seq_blocks = s[end_of_name:].split()
                sequence = ''.join(seq_blocks)
                input_sequence = ">{}\n{}".format(sequence_name, sequence)
                fasta_list.append(input_sequence)
        return fasta_list
        
    def get_column_headers_for_net_method(self, method):
        if method == 'netchop':
            titles = ('#','amino_acid','prediction_score', 'sequence_id',)
        elif method == 'netctl':
            titles = ('#','peptide','predicted_mhc_binding_affinity','rescale_binding_affinity','c_terminal_cleavage_affinity','tap_transport_efficiency','predictions_score')
        else:
            titles = ('#','peptide','mhc_prediction','tap_prediction_score','cleavage_prediction_score','combined_prediction_score','%_rank')
        return titles
    
    def remove_result_attr(self, result, method):
        """Remove all attributes from the output"""
        clean_list = []
        for tup in result:
            list_tup = list(tup)
            if method == 'netchop':
                list_tup = [ele for ele in list_tup if ele != 'Yes']
                clean_list.append(list_tup[:])
            elif method == 'netctlpan':
                removeset = set([1,2])
                clean_list.append(tuple([v for i, v in enumerate(list_tup) if i not in removeset]))
            else: clean_list.append(tup[:])
        return clean_list
    
    def print_results(self, scores, method):
        """Print scores dictionary"""
        for key, value in scores.items():
            print("{}:".format(key))
            output_col_headers = self.get_column_headers_for_net_method(method)
            print("{}".format("\t".join(output_col_headers)))
            
            clean_scores = self.remove_result_attr(value, method)
            for row in clean_scores:
                print("{}".format("\t".join(row)))

    def create_png(self, results, method, threshold, plot_name='plot'):
        """ """
        import matplotlib.pyplot as plt
        import numpy as np
        from collections import OrderedDict
        
        if type(results) == OrderedDict and len(results) == 1:
            key, values = results.items()[0]
        elif type(results) == OrderedDict and len(results) > 1:
            raise ValueError("NetChop Standalone accept one protein sequence per time. Please check your input.")
        else:
            logging.debug("Results is empty or wrong type.")
            raise ValueError("Please check your input.")

        x=[]
        y=[]
        for tup in values:
            x.append(int(tup[0]))
            if method == 'netchop':
                y.append(float(tup[2]))
            elif method == 'netctl':
                y.append(float(tup[-1]))
            else:
                y.append(float(tup[-2]))
         
        x = np.array(x)
        y = np.array(y)
        
        plt.plot(x, y)
        plt.xlabel("Position", fontsize=11)
        plt.ylabel("Score", fontsize=11)
         
        if method == "netctlpan":
            threshold = self.get_netctlpan_cutoff(values, threshold)
        plt.axhline(y=threshold, color='r', linewidth=1.5, label="Threshold")
        plt.legend(loc='upper right', prop={'size':9})
        
        plt.fill_between(x, y, np.float(threshold), where=y>np.float(threshold), color='#FFFF00', interpolate=True)
        plt.fill_between(x, y, np.float(threshold), where=y<np.float(threshold), color='#00CC00', interpolate=True)
         
        plt.grid(True)
        plt.savefig(plot_name)

    def get_netctlpan_cutoff(self, values, threshold):
        cutoff_list =[]
        score_list = []
        for value in values:
            prediction_score = float(value[-2])
            rank = float(value[-1])
            if rank <= float(threshold):
                cutoff_list.append(prediction_score)
                cutoff_list = sorted(cutoff_list)
            else: 
                score_list.append(prediction_score)
                score_list = sorted(score_list)
        
        if cutoff_list:
            cutoff = cutoff_list[0]
        else:
            indx = len(score_list) - 1
            cutoff = score_list[indx]
            cutoff += 0.001
        return cutoff
    
    def commandline_help(self):
        print """
* All available method names:
-----------------------------
netchop, netctl, netctlpan

1) Make predictions given a file containing a sequence:
-------------------------------------------------------
python predict.py -m [method-name] [input-file]
Example: python predict.py -m netchop example/test.fsa

2) Make predictions with options:
-----------------------------------------
python predict.py -m [method-name] --threshold [threshold-value] (default=0.5)
Example: python predict.py -m netchop --threshold 0.7 example/test.fsa

You can also use help option (-h or --help) for more information:
*****************************************************************
python predict.py --help
"""

    def main(self):
        
        parser = OptionParser(usage="usage: %prog -m <method_name> [options] input", version="%prog 1.0")
        
        parser.add_option("-m", "--method",
                          type="choice",
                          dest="method",
                          choices=["netchop","netctl","netctlpan"],
                          help="Select a method from available method options.",)
        parser.add_option("-t", "--threshold", 
                          dest="threshold",
                          help="Threshold value(default: netchop=0.5).")
        
        # netctl and netctlpan specific option
        parser.add_option("-w", "--tap_weight", 
                          dest="tap_weight",
                          help="NetCTL/NetCTLpan specific: weight on TAP transport efficiency (default=0.050000/0.025).")
        parser.add_option("-c", "--cleavage_weight", 
                          dest="cleavage_weight",
                          help="NetCTL/NetCTLpan specific: weight on C terminal cleavage (default: netctl=0.100000, netctlpan=0.225).")
        
        # netchop specific option
        parser.add_option("-p", "--pred_method", 
                          dest="pred_method",
                          help="Prediction method for NetChop (default=0).")
        
        # netctl specific option
        parser.add_option("-s", "--supertype", 
                          dest="supertype",
                          help="NetCTL specific: supertype(default=A1).")
        
        # netctlpan specific option
        parser.add_option("-a", "--allele", 
                          dest="allele",
                          help="NetCTLpan specific: allele(default=HLA-A02:01).")
        parser.add_option("-l", "--length", 
                          dest="length",
                          help="NetCTLpan specific: length(default=9).")
        parser.add_option("-e", "--epitope_threshold", 
                          dest="epitope_threshold",
                          help="NetCTLpan specific: epitope_threshold(default=1.000000).")
        
        # to plot or not to plot the prediction scores
        parser.add_option("-n", "--noplot",
                          action="store_true",
                          default=False, 
                          help="Do not generate a plot.")

        parser.add_option("-o", "--outputplot",
                          dest="plot_name",
                          help="The plot output path/name(default=plot).")
        
        (options, args) = parser.parse_args()       
        
        
        
        try:
            self.predict(options, args)
        except ValueError as e:
            self.commandline_help()
            parser.error(e)
    
    def args_num_validation(self, len_args):
        if len_args == 0: 
            raise ValueError("Input data file must be specified.")
        elif len_args == 1:
            return True            
        else: 
            raise ValueError("Incorrect number of arguments.")            
            
        
    
if __name__ == '__main__':
    PredictNetChop().main()

