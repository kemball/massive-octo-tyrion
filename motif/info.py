
charset = 'ACGT'
import random
def col_info(column):
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
