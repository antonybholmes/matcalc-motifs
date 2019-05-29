import sys
import numpy as np
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('file')
parser.add_argument('-t', '--threshold', default=0.5, type=float)
args = parser.parse_args()

file = args.file #'jaspar_vert_2018.pwt'


f = open(file, 'r')

for line in f:
    tokens = line.strip()[1:].split('\t')
    id = tokens[0]
    name = tokens[1]
    
    a = next(f).strip().split('\t')
    c = next(f).strip().split('\t')
    g = next(f).strip().split('\t')
    t = next(f).strip().split('\t')
    
    n = len(a)
    
    bases = [np.array([float(a[i]), float(c[i]), float(g[i]), float(t[i])]) for i in range(0, n)]
    
    score = 0
    
    for b in bases:
        print(b, file=sys.stderr)
        
        m = b.max()
        
        
        
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
        
