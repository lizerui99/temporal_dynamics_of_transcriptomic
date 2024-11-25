import sys

in_put = sys.argv[1]
out_put = sys.argv[2]

with open(in_put,'r') as inf:
    with open(out_put,'w') as of:
        for iline in inf:
            line = iline.strip().split('\t')
            if iline[0] == 'M':
                continue
            else:
                of.write(iline)