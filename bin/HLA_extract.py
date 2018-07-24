import sys
HLA_result=sys.argv[1]

HLA_str_out=sys.argv[2]

with open(HLA_result) as f:
    data=f.read()

HLA_prediction = ['HLA-'+i.strip().split(',')[0].replace('*','')[0:6] for i in data.strip().split('\n') if i.startswith('\t\t')]


HLA_str=','.join(HLA_prediction)

HLA_out_f=open(HLA_str_out,'w')

HLA_out_f.write(HLA_str)

