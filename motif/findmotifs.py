#!/usr/bin/python
from motifs import read_fasta

charset = 'ACGT'

def cheapofind(length,sequences):
    lmersets = []
    for seq in sequences:
        lmers = [seq[x:length+x] for x in range(0,len(seq)-length+1)]
        lmersets.append(set(lmers))
    return reduce(lambda x,y: x&y, lmersets)

def gen_pat(sequences,offsets,mlength):
    szip = zip(sequences,offsets)
    for x in range(0,mlength):
        letters = map(lambda z: z[0][z[1]+x],szip)
        col = [0]*len(charset)
        for num,char in enumerate(charset):
            col[num] = string.count(letters,char)







if __name__=="__main__":
    print "I should be called with some files, then find some motifs!"
