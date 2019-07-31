#!/usr/bin/env python
import shutil
import subprocess
import collections
import sys
import itertools
from optparse import OptionParser
from Bio import AlignIO
import numpy
import multiprocessing as mpr
import tempfile
from uuid import uuid4
import os

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

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

def get_ancestor(otu):
    theirfastaname = "/dev/shm/%s" % str(uuid4())
    theirtreename = "/dev/shm/%s" % str(uuid4())
    theirttdir = "/dev/shm/%s" % str(uuid4())
    theirlog = "/dev/shm/%s" % str(uuid4())
    theirrootedtreename = "/dev/shm/%s" % str(uuid4())
    theirfasta = open(theirfastaname, "w")

    cmd = ['nw_prune', '-v', tree_fp]

    for o in otu:
        theirfasta.write(">%s\n"%o)
        theirfasta.write("%s\n"%mp[o])
        cmd.append(o)
    theirfasta.close()

    with open(theirtreename, 'w') as out:
        output = subprocess.call(cmd, stdout=out)
    
    cmd = ['python', '/home/balaban/MinVar-Rooting/FastRoot.py', '-i', theirtreename, '-m', 'MV', '-o', theirrootedtreename]
    with open(theirlog, 'w') as out:
        output = subprocess.call(cmd, stdout=out)

    cmd = ['treetime', 'ancestral', '--tree', theirrootedtreename, '--aln', theirfastaname, '--outdir', theirttdir, '--gtr', 'TN93']
    with open(theirlog, 'w') as out:
        output = subprocess.call(cmd, stdout=out)
    
    rec_f = open("%s/ancestral_sequences.fasta"%theirttdir)
    for name,seq,qual in readfq(rec_f):
        if name == "NODE_0000000":
            result = seq
            break
    
    rec_f.close()

    #with open("%s/ancestral_sequences.fasta" 
    #print("theirfastaname %s" %theirfastaname)
    #print("theirrootedtreename %s" %theirrootedtreename)
    #print("theirlog %s" %theirlog)
    #print("theirttdir %s" %theirttdir)
    #print("theirttdir %s" %theirttdir)

    os.remove(theirfastaname)
    os.remove(theirrootedtreename)
    os.remove(theirtreename)
    shutil.rmtree(theirttdir)
    os.remove(theirlog)

    return seq

def collect_stats(otu):
    if len(otu) > 1:

        ancestor = get_ancestor(otu)

        dm = numpy.zeros(len(otu))
        for i in range(len(otu)):
                dm[i] = ham(mp[otu[i]], ancestor)
        mean_to_cent = numpy.mean(dm)
        max_to_cent = numpy.max(dm)
        #print(mean_to_cent,max_to_cent,len(otu))
        return (mean_to_cent,max_to_cent,len(otu))
    else:
        return (0.0,0.0,1)

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input_fp",
                      help="path to the input OTU table in gg format", metavar="FILE")
    parser.add_option("-a", "--alignment", dest="alignment_fp",
                      help="path to the input gg alignment file", metavar="FILE")
    parser.add_option("-t", "--tree", dest="tree_fp",
                      help="path to the input gg tree file", metavar="FILE")
    parser.add_option("-l", "--log", dest="log_fp",
                      help="path to log file", metavar="FILE")
    parser.add_option("-T", "--threads", dest="num_thread", default="0",
                      help="number of cores used in 0 to use all cores in the running machine", metavar="NUMBER")
# parser.add_option("-o", "--output", dest="output_fp",
   #                   help="path to the output greengenes otu table file", metavar="FILE")

    (options, args) = parser.parse_args()
    input_fp = options.input_fp
    alignment_fp = options.alignment_fp
    tree_fp = options.tree_fp
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


    #print("Reading input done")
    if num_thread:
        pool = mpr.Pool(num_thread)
    else:
        pool = mpr.Pool(mpr.cpu_count())


    #print("pool created")
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

    print("ancestral_max_max\t%.6f" % max_max) 
    print("ancestral_max_mean\t%.6f" % max_mean) 
    print("ancestral_mean_max\t%.6f" % mean_max) 
    print("ancestral_mean_mean\t%.6f" % mean_mean) 

    flog.close()

