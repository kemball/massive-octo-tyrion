#!/usr/bin/python
import random
charset = "ACGT"
def col_info(column):
    #notable: assumes bases are uniform
#normalize column for information purposes
#returns information in bits. idk whether nats are better?
    from math import log
    totes = sum (column) +0.0#lawl float
    ncol = map(lambda x: x/totes,column)
    return sum(map(lambda x:x and x*log(len(column)*x,2),ncol))

def sample_motif(motif,charset=charset):
    return "".join([charset[x] for x in map(sample_column,motif)])

def sample_column(column,charset=charset):
    total = random.randint(1,sum(column))
    for (letter,thing) in enumerate(column):
        if total<=thing:
            return letter
        else:
            total -= thing


def mot_info(motif):
    return sum(map(col_info,motif))

def create_motif(icpc, length, charset=charset):
#should return a PWM with total information icpc*length
    if icpc is not 2:
        gencol = lambda : [0.0 +random.randint(1,5) for y in range(0,len(charset))]
        motif = [gencol()  for x in range(0,length)]
        goalinfo = icpc*length
        while mot_info(motif) < goalinfo:
            if col_info(motif[0])<icpc:
                motif[0].sort()
                motif[0][-1]+= sum(motif[0])+0.0
                motif[0] = map(lambda x: x*.75,motif[0])
                random.shuffle(motif[0])
            random.shuffle(motif)
        motif = [ map(lambda x: int(x*10),col) for col in motif]
    else:
        motif = [[1]+[0]*(len(charset)-1) for x in range(0,length)]
        map(random.shuffle,motif)
        random.shuffle(motif)
    return motif

def create_benchmark(icpc=2,ml=8,sl=500,sc=10):
    sequences = [generaterandom(sl) for x in range(0,sc)]
    motif = create_motif(icpc,ml)
    offsets = []
    for x in range(0,sc):
        offset = random.randint(0,sl-ml)
        offsets.append(offset)
        q = sample_motif(motif)
        sequences[x] = sequences[x][:offset]+q+sequences[x][offset+ml:]
        if icpc == 2:
            assert(q in sequences[x])
    import os
    fileprefix = 'icpc'+str(icpc)+'ml'+str(ml)+'sl'+str(sl)+'sc'+str(sc)+'a'
    while fileprefix in os.listdir('datasets'):
        fileprefix = fileprefix[:-1]+str(chr(ord(fileprefix[-1])+1))
    os.mkdir('datasets/'+fileprefix)

    with open('datasets/' +fileprefix+'/sites.txt','w') as sitefile:
        for site in offsets:
            sitefile.write(str(site)+'\n')
    with open('datasets/'+fileprefix+'/motiflength.txt','w') as lenfile:
        lenfile.write(str(ml))
    with open('datasets/'+fileprefix+'/motif.txt','w') as mfile:
        mfile.write('>MOTIF1\t' + str(ml)+'\n')
        for col in motif:
            mfile.write("\t".join([str(x) for x in col])+'\n')
        mfile.write('<')
    with open('datasets/'+fileprefix+'/sequences.fa','w') as sfile:
        for (num,seq) in enumerate(sequences):
            sfile.write('>SEQ'+str(num)+'\n')
            sfile.write(seq+'\n')


def write_fasta(filename,strings,identifier="sequence"):
    f = open(filename,'w')
    for (num,bases) in enumerate(strings):
        f.write('>'+str(identifier) + str(num)+'\n')
        while len(bases)>=80:
            f.write(bases[0:80])
            f.write('\n')
            bases = bases[80:]
        f.write(bases)
    f.close()

def generaterandom(length=10,charset=charset):
    return "".join([ random.choice(charset) for x in range(0,length)])





if __name__=="__main__":
    for x in range(0,10):
        for icpc in [1,1.5]:
            create_benchmark(icpc=icpc)
        for ml in [6,7]:
            create_benchmark(ml=ml)
        for sc in [5,20]:
            create_benchmark(sc=sc)
        create_benchmark()
