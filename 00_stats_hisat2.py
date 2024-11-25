import os
import sys

in_put  = sys.argv[1]
out_put = sys.argv[2]
hisat2_log = os.listdir(in_put)
dict_log = {}
table_log = []

of = open(out_put,'w')
for file_name in hisat2_log:
    if file_name.split('.')[1] == 'log':
        sample_name = file_name.split('.')[0]
        with open(file_name) as in_f:
            n = 0
            for iline in in_f:
                iline = iline.strip()
                if n == 0:
                    pair_reads = iline.split(' ')[0]
                    table_log.append(pair_reads)
                    n += 1
                elif n == 3:
                    unique_reads = iline.split(')')[0]
                    unique_reads = unique_reads+")"
                    table_log.append(unique_reads)
                    n += 1
                elif n == 14:
                    mp_rate = iline.split(' ')[0]
                    table_log.append(mp_rate)
                    dict_log[sample_name] = table_log
                    print(table_log)
                    table_log = []
                    n = 0
                else:
                    n += 1

of.write('sample\tpair reads\tunique mapped(%)\tmapped rate(%)')
for sample,stats in dict_log.items():
    of.write(f'{sample}\t{stats[0]}\t{stats[1]}\t{stats[2]}\n')

of.close()
