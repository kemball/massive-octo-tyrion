import os






if __name__=="__main__":
    resfile = open('resfile.txt','w')
    resfile.write("\t".join(['icpc','ml','sl','sc','time','score'])+'\n')
    import re
    for outfile in os.listdir('datasets'):
        outfile = 'datasets/'+outfile
        stuff = filter(lambda x: x ,re.split('[a-z]',outfile[1:]))
        resfile.write("\t".join(stuff[1:]))
        fsites = []
        rsites = []
        with open(outfile+'/nsites.txt') as nsf:
            nsf.readline()
            firsttime = float(nsf.readline())
            secondtime = float(nsf.readline())
            totalt = secondtime-firsttime
            resfile.write('\t'+str(totalt)+'\t')
            fsites = nsf.readline()
        with open(outfile+'/sites.txt') as rsf:
            rsites = rsf.read()
        fsites = filter(lambda x: x,re.split('\D+',fsites))
        fsites = map(int,fsites)
        rsites = filter(lambda x: x,re.split('\D+',rsites))
        rsites = map(int,rsites)
        scores = map(lambda x: min(abs(x[0]-x[1]),8),zip(fsites,rsites))
        resfile.write(str(sum(scores))+'\t')
        resfile.write('\n')






