import pandas as pd
import numpy as np
from skbio import TabularMSA, DNA, DistanceMatrix, Protein
from skbio.sequence.distance import hamming
import os
import re
from skbio.sequence import Sequence
from numpy import isnan
from sys import stdin,stdout
import subprocess
import treeswift
import collections
from dendropy import Tree,TreeList
from Tree_extend import MPR_Tree,MV00_Tree,OGR_Tree
from Bio import Phylo

def h_distance(seq1, seq2):
    '''Modified hamming distance to include only non-gap sites'''
    myseq1 = str(seq1)
    myseq2 = str(seq2)

    degapped1 = []
    degapped2 = []

    for i in range(len(myseq1)):
        if myseq1[i] != "-" or myseq2[i] != "-":
                degapped1.append(myseq1[i])
                degapped2.append(myseq2[i])
    degapped1 = "".join(degapped1)
    degapped2 = "".join(degapped2)

    hamming_dist = hamming(Sequence(degapped1), Sequence(degapped2))
    if isnan(hamming_dist):
        return 0.0
    else:
        return hamming_dist

def get_consensus(matrix):
    output = ''
    for i in range(matrix.shape[1]):
        site = collections.Counter(matrix[:, i]).most_common(1)[0][0]
        output += site
    return output

id_seq_dict = {}
with open("/Users/jiaxingfan/Desktop/gg_13_8_otus/rep_set_aligned/99_otus.fasta") as f:
    lines = f.readlines()
    i = 0
    while i < len(lines):
        id_seq_dict[lines[i][1:-1]] = lines[i+1]
        i += 2

treeclu_th = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.06,\
               0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15]
treeclu_hamming_max = []
treeclu_hamming_mean = []
treeclu_tree_max = []
treeclu_tree_mean = []
max_max_list = []
mean_max_list = []

for th in treeclu_th:
    print(th)
    max_max_hamming = 0
    mean_max_hamming = []
    max_max_tree = 0
    mean_max_tree = []
    max_max = 0
    mean_max = []

    f = open('/Users/jiaxingfan/Desktop/run5/max/%.3f/clusters_formatted.txt' %th, 'r')
    for line in f:
        cluster = [elt.strip() for elt in line.split('\t')][1:]
        cluster_len = len(cluster)
        if cluster_len > 1:
            cmd = ['nw_prune', '-v', '/Users/jiaxingfan/Desktop/gg_13_8_otus/trees/99_otus_unannotated.tree']
            fasta_file = open('cluster_fasta.fasta', 'w')
            seq_list = []
            cmd.append(cluster[0])
            seq = id_seq_dict.get(cluster[0])
            seq_list.append(seq)
            fasta_file.write(">" + cluster[0] + "\n")
            fasta_file.write(seq + "\n")
            # consensus start
            mat = np.array(list(id_seq_dict.get(cluster[0]))[:-1])
            for j in range(1, len(cluster)):
                node = cluster[j]
                cmd.append(node)
                seq = id_seq_dict.get(node)
                mat = np.vstack((mat, np.array(list(seq[:-1]))))
                seq_list.append(seq)
                fasta_file.write(">" + node + "\n")
                fasta_file.write(seq + "\n")
            consensus = get_consensus(mat)
            temp_max = 0
            fasta_file.close()
            for j in range(len(cluster)):
                temp_dist = h_distance(consensus, ''.join(mat[j, :]))
                temp_max = max(temp_dist, temp_max)
            max_max = max(max_max, temp_max)
            mean_max.append(temp_max)
            # consensus end

            # ancestral reconstruction start
            # check identical sequences
            if len(seq_list) != len(set(seq_list)):
                print("identical sequences!")
                print(seq_list)
                continue

            temp_max_hamming = 0
            temp_max_tree = 0
            with open('cluster.tree', 'w') as out:
                output = subprocess.call(cmd, stdout=out)

            # apply re-rooting
            tree_file = open('cluster.tree', 'r')
            for tree_line in tree_file:
                tree = Tree.get(data = tree_line, schema = "newick", preserve_underscores = True)
                # min variance
                a_tree = MV00_Tree(ddpTree = tree)
                # midpoint
                #a_tree = MPR_Tree(ddpTree = tree)
            tree_file.close()

            re_rooted_tree = open('re_rooted_tree.nwk', 'w')
            re_rooted_tree.write(a_tree.ddpTree.as_string("newick"))
            re_rooted_tree.close()

            # apply ancestral reconstruction
            os.system('treetime ancestral --tree re_rooted_tree.nwk --aa --aln cluster_fasta.fasta --outdir ancestral_results >useless_output.txt')
            ances = open('ancestral_results/ancestral_sequences.fasta', 'r')
            lines = ances.readlines()
            centroid = lines[1][:-1]
            cent_label = lines[0][1:-1]
            ####### need revise #######################################################
            for line_num in range(2, 130):
                centroid += lines[line_num][:-1]
            ances.close()

            for member in seq_list:
                temp_dist_hamming = h_distance(centroid, member)
                temp_max_hamming = max(temp_dist_hamming, temp_max_hamming)
            max_max_hamming = max(max_max_hamming, temp_max_hamming)
            mean_max_hamming.append(temp_max_hamming)

            # check if there are identical sequences
            # get tree distance
            #Phylo.convert('ancestral_results/annotated_tree.nexus', 'nexus', 'ancestral_results/annotated_tree.newick', 'newick')
            #clu_tree = treeswift.read_tree_newick('ancestral_results/annotated_tree.newick')
            with open("ancestral_results/annotated_tree.nexus") as nexus_file:
                nexus_line = nexus_file.readline()
                nexus_line = nexus_file.readline()
                nexus_line = nexus_file.readline()
                nexus_line = nexus_file.readline()
                nexus_line = nexus_file.readline()
                nexus_line = nexus_file.readline()
                nexus_line = nexus_file.readline()
            nexus_line = re.sub("[\[].*?[\]]", "", nexus_line)[12:]
            clu_tree = treeswift.read_tree(nexus_line, "Newick")
            #lab_to_node = clu_tree.label_to_node()
            #cent_node_node = lab_to_node.get(cent_label)
            #print(cent_node_node)
            #print(lab_to_node)
            #print(cent_label)
            dis_tuples = clu_tree.distances_from_root()
            for items in dis_tuples:
                cur_dist = items[1]
                temp_max_tree = max(temp_max_tree, cur_dist)
            #for i in range(cluster_len):
            #    cur_dist = clu_tree.distance_between(cent_node_node, lab_to_node.get(cluster[i]))
            #    temp_max_tree = max(temp_max_tree, cur_dist)

            mean_max_tree.append(temp_max_tree)
            max_max_tree = max(max_max_tree, temp_max_tree)
            # ancestral reconstruction end

    # consensus
    max_max_list.append(max_max)
    mean_max_list.append(np.mean(mean_max))
    # ancestral reconstruction
    treeclu_hamming_max.append(max_max_hamming)
    treeclu_hamming_mean.append(np.mean(mean_max_hamming))
    treeclu_tree_max.append(max_max_tree)
    treeclu_tree_mean.append(np.mean(mean_max_tree))
    f.close()
    with open('/Users/jiaxingfan/Desktop/gg_max_out.txt', 'a') as out_f:
        f3.write(' '.join(list(map(str, treeclu_hamming_max))))
        f3.write(' '.join(list(map(str, treeclu_hamming_mean))))
        f3.write(' '.join(list(map(str, treeclu_tree_max))))
        f3.write(' '.join(list(map(str, treeclu_tree_mean))))

print('treeclu_hamming_max = ', treeclu_hamming_max)
print('treeclu_hamming_mean = ', treeclu_hamming_mean)
print('treeclu_tree_max = ', treeclu_tree_max)
print('treeclu_tree_mean = ', treeclu_tree_mean)

