#!/usr/bin/env python
import sys
import itertools
import treeswift
from optparse import OptionParser


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input_fp",
                      help="path to the input OTU table to be extended in gg format", metavar="FILE")
    parser.add_option("-p", "--previous", dest="previous_fp",
                      help="path to OTU table higher in the hierarchy", metavar="FILE")
    parser.add_option("-o", "--output", dest="output_fp",
                      help="path to the expanded OTU table file", metavar="FILE")

    (options, args) = parser.parse_args()
    input_fp = options.input_fp
    output_fp = options.output_fp
    previous_fp = options.previous_fp

    #print("Reading input")
    fin = open(input_fp, "r")
    otus = map(lambda x: x.strip().split('\t')[1:], fin.readlines())
    fin.close()
    fp = open(previous_fp, "r")
    otus_p = map(lambda x: x.strip().split('\t')[1:], fp.readlines())
    fp.close()
    rep_seqs= dict([(otu[0], otu[1:]) for otu in otus_p])
    new_otus = [list(itertools.chain.from_iterable([[i] + rep_seqs[i] for i in otu])) for otu in otus]
    fout = open(output_fp, "w")
    count = 0
    for i in new_otus:
        fout.write(str(count) + '\t' + '\t'.join(i) + '\n')
        count +=1
    fout.close()

