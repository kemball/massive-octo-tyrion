#!/usr/bin/python
import random
charset = "ACGT"
def information(column):
    #notable: assumes bases are uniform
#normalize column for information purposes
#returns information in bits. idk whether nats are better?
    from math import log
    totes = sum (column)
    ncol = map(lambda x: x/totes,column)
    return sum(map(lambda x:x and x*log(len(charset)*x,2),ncol))

def create_motif(icpc, length, charset=charset):
#should return a PWM with total information icpc*length
    gencol = lambda : [random.randint(1,10) for y in range(0,len(charset))]
    motif = [ gencol() for x in range(0,length)]
    goalinfo = icpc*length
    random.shuffle(motif)
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
