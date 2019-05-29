import sys
import numpy as np
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('file')
parser.add_argument('-t', '--threshold', default=0.5, type=float)
args = parser.parse_args()

file = args.file #'jaspar_vert_2018.pwt'


f = open(file, 'r')

f.readline()
f.readline()

for line in f:
    tokens = line.strip().split('\t')
    
    id = tokens[0]
    name = tokens[1]
    bases = [np.array([float(x) for x in bases.split(',')]) for bases in tokens[3].split(';')]
    
    for i in range(0, len(bases)):
        bases[i] = (bases[i] / bases[i].sum()).round(3)
    
    score = 0
    
    for b in bases:
        m = b.max()
        
        #print(m, file=sys.stderr)
        
        l = np.log2(m / 0.25)
        
        score += l
    
    #print(score, file=sys.stderr)
    #print(bases, file=sys.stderr)
    
    threshold_score = score * args.threshold
    
    print('>{}\t{}\t{}'.format(id, name, threshold_score))
    
    for b in bases:
        print('\t'.join([str(s) for s in b]))
    
    #break
    
f.close()
        
