import re
import os
from subprocess import Popen, PIPE
from collections import OrderedDict
from curses.ascii import isdigit

TMP_DIR_PATH = '/tmp'
EXECUTABLE_NAME = 'netChop'
EXECUTABLE_DIR_PATH = os.path.abspath(os.path.dirname(__file__))
EXECUTABLE_FULL_PATH = os.path.join(EXECUTABLE_DIR_PATH, 'bin', EXECUTABLE_NAME)

# set the environment variables to the 'netchop-3.1' and temporary directory
os.environ["NETCHOP"] = EXECUTABLE_DIR_PATH
os.environ["TMPDIR"] = TMP_DIR_PATH

def predict(fasta_list, threshold='0.500000', method='0', short=False):
    """ @brief: runs executable, parses result and returns a score dictionary
        @params: list of input fasta sequences, threshold(default=0.500000), method(default=0) and short(default=False)
    """
    import itertools
    import operator
    import tempfile
    
    # check if input is a list
    assert type(fasta_list) is list, "Input file must be a list of fasta sequence(s)."
        
    # write a temporary file from a fasta_list items
    tmpfile = tempfile.NamedTemporaryFile()
    for fasta_sequence in fasta_list:
        tmpfile.write("{}\n".format(fasta_sequence))
    tmpfile.seek(0)
    
    cmd = [EXECUTABLE_FULL_PATH, '-t', threshold, '-v', method, tmpfile.name]
    if short: cmd.append('-s')
    process = Popen(cmd, stdout=PIPE)
    (stdout_value, stderr_value) = process.communicate()
    
    scores = OrderedDict()
    if short is False:
        output = stdout_value.splitlines()
        
        # list_res is list of tuples of prediction scores without the comment lines and column headers
        list_res = [res.split() for res in output if not res.startswith('#') and res.split() and res.split()[0].isdigit()]
        
        # Remove 'Asigned Prediction' ('S' for prediction > threshold, '.' otherwise) attributes from the output
        list_res = remove_result_attr(list_res)
    
        # group the list by the last tuple element
        for key, group in itertools.groupby(list_res, operator.itemgetter(-1)):
            filtered_group_list = []
            for lis in list(group):
                filtered_group_list.append(lis)
            scores[key] = filtered_group_list
    else:
        # leave the first part of the prediction result
        seqs = re.split('NetChop\s.+', stdout_value)[1]
        seqs = re.split('Number of cleavage sites.+Number of amino acids.+', seqs)
        seqs = [re.sub('\-+', '', res).strip().split() for res in seqs]
        for seq in seqs:
            if seq: scores[seq[1]] = seq[2:]
    
    return scores

def remove_result_attr(result):
    """Remove all attributes from the output"""
    filtered_list = []
    for tup in result:
        lis = list(tup)
        removeset = set([2])
        filtered_list.append(tuple([v for i, v in enumerate(lis) if i not in removeset]))
    return filtered_list
