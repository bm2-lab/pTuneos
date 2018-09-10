import os
from subprocess import Popen, PIPE
from collections import OrderedDict

TMP_DIR_PATH = '/tmp'
EXECUTABLE_NAME = 'netCTLpan'
EXECUTABLE_DIR_PATH = os.path.relpath(os.path.dirname(__file__))
EXECUTABLE_FULL_PATH = os.path.join(EXECUTABLE_DIR_PATH, EXECUTABLE_NAME)

# set the environment variables to the 'netCTLpan-1.1' and temporary directory
os.environ["NETCTLPAN"] = EXECUTABLE_DIR_PATH
os.environ["TMPDIR"] = TMP_DIR_PATH

def predict(fasta_list, allele='HLA-A02:01', length='9', tap_weight='0.025000', cleavage_weight='0.225000', epitope_threshold='1.000000', threshold='-99.900002', sort='-1'):
    """ @brief: runs executable, parses result and returns a score dictionary
        @params: list of input fasta sequences, allele(default='HLA-A02:01'), length(default='9'), tap_weight(default='0.025000'), cleavage_weight(default='0.225000'), 
                 epitope_threshold(default='1.000000'), threshold(default='-99.900002') and sort(default='-1')
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
    
    cmd = [EXECUTABLE_FULL_PATH, '-a', allele, '-l', length, '-wt', tap_weight, '-wc', cleavage_weight, '-ethr', epitope_threshold, '-thr', threshold, '-s', sort, tmpfile.name]
    process = Popen(cmd, stdout=PIPE)
    (stdout_value, stderr_value) = process.communicate()
     
    scores = OrderedDict()
    output = stdout_value.splitlines()
     
    # list_res is list of tuples of prediction scores without the comment lines and column headers
    list_res = [tuple(res.split()) for res in output if not res.startswith('#') and res.split() and res.split()[0].isdigit()]
    
    # remove any non-score elements from the tuple
    list_res = cleanup(list_res)
    
    # group the list by the second tuple element
    for key, group in itertools.groupby(list_res, operator.itemgetter(1)):
        scores[key] = list(group)
        
    return scores

def cleanup(result):
    """ remove anything besides the 8 elements/scores from the output and update numbering to start from 1 """
    cleanded_list = []
    for tup in result:
        lis = list(tup)
        lis[0] = str(int(lis[0]) + 1)
        cleanded_list.append(tuple([i for i in lis[:9]]))
    return cleanded_list
    