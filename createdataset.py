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
def mot_info(motif):
    return sum(map(col_info,motif))

def create_motif(icpc, length, charset=charset):
#should return a PWM with total information icpc*length
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
    return motif

def create_benchmark(icpc=2,ml=8,sl=500,sc=10):
    charset = charset
    sequences = [generaterandom(length=sl) for x in range(0,sc)]
    pass

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
    print "I create datasets in datasets/something!"
