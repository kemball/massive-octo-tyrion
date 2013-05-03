#!/usr/bin/python
from motifs import read_fasta

def cheapofind(length,sequences):
    lmersets = []
    for seq in sequences:
        lmers = [seq[x:length+x] for x in range(0,len(seq)-length+1)]
        lmersets.append(set(lmers))
    return lmersets
    return reduce(lambda x,y: x&y, lmersets)






if __name__=="__main__":
    print "I should be called with some files, then find some motifs!"
