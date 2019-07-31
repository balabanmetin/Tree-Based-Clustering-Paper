#!/usr/bin/env python
import collections
import sys
import itertools
from optparse import OptionParser
from Bio import AlignIO
import numpy
import multiprocessing as mpr


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

def collect_stats(otu):
    if len(otu) > 1:
        otu_seqs = []
        for sq in otu:
            otu_seqs.append(mp[sq])
        k = numpy.array(list(map(lambda x: list(x), otu_seqs)))
        consensus = get_consensus(k)

        dm = numpy.zeros(len(otu))
        for i in range(len(otu)):
                dm[i] = ham(mp[otu[i]], consensus)
        mean_to_cent = numpy.mean(dm)
        max_to_cent = numpy.max(dm)
        return (mean_to_cent,max_to_cent,len(otu))
    else:
        return (0.0,0.0,1)

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input_fp",
                      help="path to the input OTU table in gg format", metavar="FILE")
    parser.add_option("-a", "--alignment", dest="alignment_fp",
                      help="path to the input gg alignment file", metavar="FILE")
    parser.add_option("-l", "--log", dest="log_fp",
                      help="path to log file", metavar="FILE")
    parser.add_option("-T", "--threads", dest="num_thread", default="0",
                      help="number of cores used in 0 to use all cores in the running machine", metavar="NUMBER")
# parser.add_option("-o", "--output", dest="output_fp",
   #                   help="path to the output greengenes otu table file", metavar="FILE")

    (options, args) = parser.parse_args()
    input_fp = options.input_fp
    alignment_fp = options.alignment_fp
    log_fp = options.log_fp
    num_thread = int(options.num_thread)


    #print("Reading input")
    fin = open(input_fp, "r")
    otus = map(lambda x: x.strip().split('\t')[1:], fin.readlines())
    fin.close()
    alignment = AlignIO.read(open(alignment_fp), 'fasta')
    mp = dict([(record.id, record.seq.__str__()) for record in alignment])

    flog = open(log_fp, "w")

    weights = 0
    max_max = 0
    max_mean = 0
    mean_max = 0
    mean_mean = 0


    if num_thread:
        pool = mpr.Pool(num_thread)
    else:
        pool = mpr.Pool(mpr.cpu_count())


    results = pool.map(collect_stats, otus)

    for i,j,k in results:
        if k > 1:
            weights += 1
            max_max = max(max_max, j)
            max_mean = max(max_mean, i)
            mean_max += j
            mean_mean += i
    mean_max = mean_max * 1.0 / weights
    mean_mean = mean_mean * 1.0 / weights

    print("consensus_max_max\t%.6f" % max_max) 
    print("consensus_max_mean\t%.6f" % max_mean) 
    print("consensus_mean_max\t%.6f" % mean_max) 
    print("consensus_mean_mean\t%.6f" % mean_mean) 

    flog.close()

