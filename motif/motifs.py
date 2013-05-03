#!/usr/bin/python
def read_fasta(filename):
    with open(filename,'r') as f:
        sc = 0
        seq=[""]
        for line in f:
            if '>' in line:
                sc = sc+1
                seq.append("")
                continue
            if sc:
                seq[sc] = seq[sc]+line.rstrip()
        return seq[1:]







if __name__=="__main__":
    print "I should be called with some files, then find some motifs!"
