#!/usr/bin/python
from motifs import read_fasta
from info import *
#debug
seq = read_fasta('../datasets/test/sequences.fa')


def cheapofind(length,sequences):
    lmersets = []
    for seq in sequences:
        lmers = [seq[x:length+x] for x in range(0,len(seq)-length+1)]
        lmersets.append(set(lmers))
    return reduce(lambda x,y: x&y, lmersets)

def gen_pat(sequences,offsets,mlength):
    szip = zip(sequences,offsets)
    motif = []
    for x in range(0,mlength):
        letters = map(lambda z: z[0][z[1]+x],szip)
        col = [0]*len(charset)
        for num,char in enumerate(charset):
            col[num] = letters.count(char)
        motif.append(col)
    return motif

def score_motif(sequences,motif,offsets):
    motifchunks = []
    for seqnum,offset in enumerate(offsets):
        motifchunks.append(sequences[seqnum][offset:offset+len(motif)])
    score = 0.0
    for chunk in motifchunks:
        score += score_chunk(chunk,motif)
    return score/float(len(motif))

def score_chunk(chunk,motif):
    nchunk = map(lambda x: charset.find(x),chunk)
    return sum([col[num]/float(sum(col)) for num,col in zip(nchunk,motif)])
        nchunk = map(lambda x: charset.find(x),chunk)
        score +=sum([col[num]/float(sum(col))for num,col in zip(nchunk,motif)])






if __name__=="__main__":
    print "I should be called with some files, then find some motifs!"
