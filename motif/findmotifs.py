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

def lfind(s,sub):
    l = []
    start = 0
    while s.find(sub,start) is not -1:
        l.append(s.find(sub,start))
        start +=s.find(sub,start)+len(sub)
    return l


def offset_gen(sequence,motif):
    minfo = mot_info(motif)
    keepnum = len(sequence)/2**minfo+1.0
    keepnum = int(keepnum)
    overlap = []
    while len(overlap)==0:
        threemers = lambda splot: [splot[x:x+3] for x in range(0,len(splot)-2)]
        sample = sample_motif(motif)
        setsamp = set(threemers(sample))
        seqthreemers = threemers(sequence)
        overlap = set(seqthreemers)&setsamp
    offsets = lfind(sequence,random.choice(list(overlap)))
    chunkbetter = lambda a,b: score_chunk(sequence[a:a+len(motif)],motif)-score_chunk(sequence[b:b+len(motif)],motif)
    scored_offsets = sorted(offsets,cmp=chunkbetter)
    return random.choice(scored_offsets[:keepnum])

def gibbs_iter(seq,length):
    off = [random.randint(0,len(s)-length) for s in seq]
#loop
    szip = zip(seq,off)
    random.shuffle(szip)
    seq,off = unzip(szip)
    m=gen_pat(seq[1:],off[1:0],length)



if __name__=="__main__":
    print "I should be called with some files, then find some motifs!"
