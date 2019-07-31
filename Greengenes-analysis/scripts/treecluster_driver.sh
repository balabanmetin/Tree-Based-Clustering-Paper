#!/bin/bash

source activate three

    mkdir $1
    python ~/TreeCluster/TreeCluster.py -i /oasis/projects/nsf/uot138/balaban/treecluster/gg_13_8_otus/trees/99_otus_unannotated.tree -o $1/clusters.txt -m $2 -t $3

    python /home/balaban/clusterproj/scripts/treecluster2gg.py -i $1/clusters.txt -o $1/clusters_formatted.txt
