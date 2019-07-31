#!/usr/bin/env python
import collections
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

def get_consensus(matrix):
    output = ''
    for i in range(matrix.shape[1]):
        site = collections.Counter(matrix[:, i]).most_common(1)[0][0]
        output += site
    return output

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input_fp",
                      help="path to the input OTU table in gg format", metavar="FILE")
    parser.add_option("-a", "--alignment", dest="alignment_fp",
                      help="path to the input gg alignment file", metavar="FILE")
    parser.add_option("-l", "--log", dest="log_fp",
                      help="path to log file", metavar="FILE")
   # parser.add_option("-o", "--output", dest="output_fp",
   #                   help="path to the output greengenes otu table file", metavar="FILE")

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
        if len(otu) > 1:
            otu_seqs = []
            for sq in otu:
                otu_seqs.append(mp[sq])
            k = numpy.array(list(map(lambda x: list(x), otu_seqs)))
            consensus = get_consensus(k)

            dm = numpy.zeros(len(otu))
            for i in range(len(otu)):
                    dm[i] = ham(mp[otu[i]], consensus)
            maxi = numpy.max(dm)
            flog.write(str(len(dm)) + '\t' + str(maxi) + '\n')
            ind += maxi
            weights += 1
            #print(ind * 1.0 / weights)
    ind = ind * 1.0 / weights
    print(ind)
    flog.close()

