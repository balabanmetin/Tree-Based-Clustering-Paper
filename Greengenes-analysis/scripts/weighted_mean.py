#!/usr/bin/env python
import sys
import itertools
import treeswift
from optparse import OptionParser


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input_fp",
                      help="path to the input OTU table in gg format", metavar="FILE")
    parser.add_option("-t", "--tree", dest="tree_fp",
                      help="path to the input gg tree file", metavar="FILE")
    parser.add_option("-l", "--log", dest="log_fp",
                      help="path to log file", metavar="FILE")

    (options, args) = parser.parse_args()
    input_fp = options.input_fp
    tree_fp = options.tree_fp
    log_fp = options.log_fp

    #print("Reading input")
    fin = open(input_fp, "r")
    otus = map(lambda x: x.strip().split('\t')[1:], fin.readlines())
    fin.close()
    flog = open(log_fp, "w")
    #print("Reading tree")
    #ftree = open(tree_fp,"r")
    tree = treeswift.read_tree_newick(tree_fp)
    #ftree.close()
    ind = 0
    weights = 0
    for otu in otus:
        if len(otu) == 1:
            weights += 1
            flog.write(str(1) + '\t' + str(0) + '\n')
        if len(otu) > 1:
            subt = tree.extract_tree_with(otu, suppress_unifurcations = True)
            dm = subt.distance_matrix()
            sum = 0
            for k1,v1 in dm.items():
                for k2, v2 in v1.items():
                    sum += v2
            flog.write(str(len(dm)) + '\t' + str(sum) + '\n')
            sum = sum * 1.0 / len(dm)
            ind += sum
            weights += len(dm)
            #print(ind * 1.0 / weights)
    ind = ind * 1.0 / weights
    print(ind)
    flog.close()


