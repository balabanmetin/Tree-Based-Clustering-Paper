#!/usr/bin/env python
import sys
import itertools
from optparse import OptionParser
from Bio import AlignIO
import numpy

def ham(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    mm=0
    sup=0
    for el1, el2 in zip(s1, s2):
        if el1=='-' and el2=='-':
            continue
        sup +=1
        if el1 != el2:
            mm +=1
    return mm*1.0/sup

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input_fp",
                      help="path to the input OTU table in gg format", metavar="FILE")
    parser.add_option("-a", "--alignment", dest="alignment_fp",
                      help="path to the input gg alignment file", metavar="FILE")
    parser.add_option("-l", "--log", dest="log_fp",
                      help="path to log file", metavar="FILE")
    #parser.add_option("-o", "--output", dest="output_fp",
    #                  help="path to the output greengenes otu table file", metavar="FILE")

    (options, args) = parser.parse_args()
    input_fp = options.input_fp
    alignment_fp = options.alignment_fp
    log_fp = options.log_fp


    #print("Reading input")
    fin = open(input_fp, "r")
    otus = map(lambda x: x.strip().split('\t')[1:], fin.readlines())
    fin.close()
    alignment = AlignIO.read(open(alignment_fp), 'fasta')
    mp = dict([(record.id, record.seq.__str__()) for record in alignment])

    flog = open(log_fp, "w")

    ind = 0
    weights = 0
    for otu in otus:
        if len(otu) >= 1:
            dm = numpy.zeros((len(otu), len(otu)))
            for i in range(len(otu)):
                for j in range(i+1, len(otu)):
                    dm[i][j] = ham(mp[otu[i]], mp[otu[j]])
                    dm[j][i] = dm[i][j]
            sum = numpy.sum(dm)
            flog.write(str(len(dm)) + '\t' + str(sum) + '\n')
            sum = sum * 1.0 / len(dm)
            ind += sum
            weights += len(dm)
            #print(ind * 1.0 / weights)
    ind = ind * 1.0 / weights
    print(ind)
    flog.close()

