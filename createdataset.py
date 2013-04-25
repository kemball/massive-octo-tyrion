#!/usr/bin/python
import random
from math import log
charset = "ACGT"
def information(column):
    #notable: assumes bases are uniform
    return sum(map(lambda x:x and x*log(len(charset)*x),column))

def create_motif(icpc, length, charset=charset):
#should return a PWM with
    pass

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
