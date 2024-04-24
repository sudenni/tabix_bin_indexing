""" Tabix binning algorithm in Python -- https://samtools.github.io/hts-specs/tabix.pdf """
import argparse
import numpy as np

def reg2bin(beg, end):
    """Caculate the bin number"""
    end -= 1
    if beg >> 14 == end >> 14:
        return ((1<<15)-1)/7 + (beg>>14)
    if beg >> 17 == end >> 17:
        return ((1<<12)-1)/7 + (beg>>17)
    if beg >> 20 == end >> 20:
        return ((1<<9)-1)/7 + (beg>>20)
    if beg >> 23 == end >> 23:
        return ((1<<6)-1)/7 + (beg>>23)
    if beg >> 26 == end >> 26:
        return ((1<<3)-1)/7 + (beg>>26)
    return 0

def reg2bins(rbeg, rend, max_bin):
    """List of bins overlapping a region"""
    # i = 0
    li = np.zeros(int(max_bin)).tolist()
    rend -= 1
    i=1
    li[i] = 0
    for k in range((1 + (rbeg>>26)), (1 + (rend>>26))):
        i+=1
        li[i] = k
    for k in range((9 + (rbeg>>23)), (9 + (rend>>23))):
        i+=1
        li[i] = k
    for k in range((73 + (rbeg>>20)), (73 + (rend>>20))):
        i+=1
        li[i] = k
    for k in range((585 + (rbeg>>17)), (585 + (rend>>17))):
        i+=1
        li[i] = k
    for k in range((4681 + (rbeg>>14)), (4681 + (rend>>14))):
        i+=1
        li[i] = k
    return i


parser = argparse.ArgumentParser(description='Binning')
parser.add_argument('begining', type = int)
parser.add_argument('end', type = int)

args = parser.parse_args()
print(reg2bin(args.begining, args.end))

max_bin = ((1<<18)-1)/7
print(reg2bins(args.begining, args.end, max_bin))