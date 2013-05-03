#!/usr/bin/python
from motifs import read_fasta
from info import *
from math import log
#debug
seq = read_fasta('../datasets/test/sequences.fa')


def cheapofind(length,sequences):
    lmersets = []
    for seq in sequences:
        lmers = [seq[x:length+x] for x in range(0,len(seq)-length+1)]
        lmersets.append(set(lmers))
    return list(reduce(lambda x,y: x&y, lmersets))

def gen_pat(sequences,offsets,mlength):
    szip = zip(sequences,offsets)
    motif = []
    for x in range(0,mlength):
        letters = []
        for z in szip:
            letters.append(z[0][z[1]+x])
        letters = map(lambda z: z[0][z[1]+x],szip)
        col = [0]*len(charset)
        for num,char in enumerate(charset):
            col[num] = letters.count(char)
        motif.append(col)
    return motif

def score_motif(sequences,motif,offsets):
    return mot_info(gen_pat(sequences,offsets,len(motif)))
    motifchunks = []
    for seqnum,offset in enumerate(offsets):
        motifchunks.append(sequences[seqnum][offset:offset+len(motif)])
    score = 0.0
    for chunk in motifchunks:
        score += score_chunk(chunk,motif)
    return score

def score_chunk(chunk,motif):
    nchunk = map(lambda x: charset.find(x),chunk)
    rnum =5000
    return reduce(lambda x,y: (1+rnum*x)*(1+rnum*y),([col[num]/float(sum(col)) for num,col in zip(nchunk,motif)]))

def lfind(s,sub):
    l = []
    start = 0
    while s.find(sub,start) is not -1:
        l.append(s.find(sub,start))
        start +=s.find(sub,start)+len(sub)
    return l


def offset_gen(sequence,motif):
    minfo = mot_info(motif)
    mers = list(set([sequence[x:x+len(motif)] for x in range(0,len(sequence)-len(motif))]))
    possibles = [(mer,score_chunk(mer,motif)) for mer in mers]
    possibles.sort(key=lambda x: x[1])
    total = max([poss[1] for poss in possibles])
    key = random.random()*total
    displace=-1
    while key >0:
        displace+=1
        key -= possibles[displace][1]
    return random.choice(lfind(sequence,possibles[displace][0]))

def bootstrap(seq,length,min = 1):
    s= set([])
    thresh = length
    while len(s) < min:
        s = cheapofind(thresh,seq)
        thresh -= 1
    s = list(s)
    random.shuffle(s)
    for el in s:
        poss = [lfind(thisseq,el) for thisseq in seq]
        poss2 = [filter(lambda small: small<(len(seq[0])-length),p) for p in poss]
        if [] in poss2:
            continue
        yield [random.choice(p) for p in poss2]


def gibbs_iter(seq,length,iters=5000):
    off = random.choice(list(bootstrap(seq,length,1)))
    nseq = seq[:]
    for zx in xrange(0,iters):
        szip = zip(nseq,off)
        random.shuffle(szip)
        nseq,off = zip(*szip)
        nseq = list(nseq)
        off = list(off)
        m=gen_pat(nseq[1:],off[1:],length)
        off[0] = offset_gen(nseq[0],m)
        if zx%(len(seq[0])/2) == 0:
            print mot_info(m)
            if mot_info(m) >= 1.5*length:
                soff = []
                for s in seq:
                    soff.append(off[nseq.index(s)])
                return soff
            smot = score_motif(seq,m,off)
            if max(off)<len(seq[0])-8 and score_motif(seq,m,[o+1 for o in off])>smot:
                    off = [o+1 for o in off]
            if min(off)>0 and score_motif(seq,m,[o-1 for o in off])>smot:
                    off = [o-1 for o in off]
            if False:
                print "mot_info(m) ("+str(mot_info(m))+") too small... restarting"
                off = [random.randint(0,len(seq[0])-length) for s in seq]
                continue
    #now I have to put things back in order? wooooork
    soff = []
    for s in seq:
        soff.append(off[nseq.index(s)])
    return soff

def chunks(seq,off,length):
    return [s[o:o+length] for s,o in zip(seq,off)]

def test():
    off = gibbs_iter(seq,8)
    for c in chunks(seq,off,8):
        print c
    print mot_info(gen_pat(seq,off,8))
    return off


if __name__=="__main__":
    print "I should be called with some files, then find some motifs!"
