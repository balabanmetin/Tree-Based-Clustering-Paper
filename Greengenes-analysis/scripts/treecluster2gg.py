#!/usr/bin/env python
import sys
import itertools
from optparse import OptionParser


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input_fp",
                      help="path to the input treecluster output file", metavar="FILE")
    parser.add_option("-o", "--output", dest="output_fp",
                      help="path to the output greengenes otu table file", metavar="FILE")

    (options, args) = parser.parse_args()
    input_fp = options.input_fp
    output_fp = options.output_fp

    fin = open(input_fp, "r")
    fout = open(output_fp, "w")
    tobewasted = fin.readline()
    lines = map(lambda x: x.strip().split('\t'), fin.readlines())
    lines_sorted = sorted(lines, key=lambda x: x[1])
    counter = 0
    for key, group in itertools.groupby(lines_sorted, lambda x: x[1]):
        if key == "-1":
            for thing in group:
                fout.write(str(counter) + '\t'+thing[0] + "\n")
                counter += 1
        else:
            things = [i[0] for i in group]
            fout.write(str(counter) + '\t' + '\t'.join(things) + '\n')
            counter += 1
            #for thing in group:
            
            #fout.write(thing[0]) 
            #fout.write("\n")
    fout.close()
    fin.close()
