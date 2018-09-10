import os
import re
from subprocess import Popen, PIPE
from collections import OrderedDict

TMP_DIR_PATH = '/tmp'
EXECUTABLE_NAME = 'netCTL'
EXECUTABLE_DIR_PATH = os.path.abspath(os.path.dirname(__file__))
# Replace a longer 'netctl_1_1' directory path inside the virtual environment to a short symlinked directory. 
# This hack is applied for a bug which is related (strangely) to character limit in the path
EXECUTABLE_DIR_PATH = re.sub('^/opt/python-virtualenvs/djangotools-website/lib/python2.7/site-packages/netctl_1_1', '/opt/netctl_1_1', EXECUTABLE_DIR_PATH)
EXECUTABLE_FULL_PATH = os.path.join(EXECUTABLE_DIR_PATH, 'bin', EXECUTABLE_NAME)

# set the environment variables to the 'netctl-1.1' and temporary directory
os.environ["NETCTL"] = EXECUTABLE_DIR_PATH
os.environ["TMPDIR"] = TMP_DIR_PATH

def predict(fasta_list, tap_weight='0.050000', mhc_weight='1.000000', cleavage_weight='0.100000', supertype='A1', sort_output='-1', threshold='0.750000'):
    """ @brief: runs executable, parses result and returns a score dictionary
        @params: list of input fasta sequences, tap_weight(default='0.050000'), mhc_weight(default='1.000000'), cleavage_weight(default='0.100000'),
                 supertype(default='A1'), sort_output(default='-1') and threshold(default='0.750000')
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
        
    cmd = [
        EXECUTABLE_FULL_PATH, '-thr', threshold, '-wc', cleavage_weight, '-wt', tap_weight, 
        '-s', supertype, '-we', mhc_weight, '-sort', sort_output, tmpfile.name
    ]
    
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    process_status_code = p.wait()
    analysis_results, ignored_stderr = p.communicate()
    tmpfile.close()

    scores = OrderedDict()
    output = analysis_results.splitlines()
    
    if process_status_code != 0:
        msg = 'Error calling netctl-1.1 executable:\n{}'.format(' '.join(cmd))
        raise Exception(msg)
    
    # list_res is list of tuples of prediction scores without the comment lines and column headers
    list_res = [res.split() for res in output if not res.startswith('#') and res.split() and res.split()[0].isdigit()]
    
    # Remove 'ID','pep','aff','aff_rescale','cle','tap','COMB' attributes from the output
    list_res = remove_result_attr(list_res)
    
    # group the list by the second tuple element
    for key, group in itertools.groupby(list_res, operator.itemgetter(1)):
        filtered_group_list = []
        list_of_list = map(list, list(group))
        for lis in list_of_list:
            # remove the last list element(input name)
            del lis[1]
            filtered_group_list.append(tuple(lis))
        scores[key] = filtered_group_list
    return scores

def remove_result_attr(result):
    """Remove all attributes from the output"""
    filtered_list = []
    for tup in result:
        lis = list(tup[:15])
        removeset = set([1,3,5,7,9,11,13])
        filtered_list.append(tuple([v for i, v in enumerate(lis) if i not in removeset]))
    return filtered_list
