#!/usr/bin/env bash
if [ "$#" -ne 4 ]; then
    echo "USAGE: $0 <tree> <mode> <threshold> <num>"
fi
tree=$1
mode=$2
thresh=$3
num=$4
clusters="$tree.clusters.txt"
trans="$(echo $tree | rev | cut -d'/' -f2- | rev)/$(echo $tree | rev | cut -d'/' -f1 | rev | cut -d'.' -f1).transmissions.txt.gz"
time8="$clusters.time8"
time9="$clusters.time9"
growth="$tree.growth.txt"
people="$tree.people.txt"
zcat $tree | python3 MinVar-Rooting/FastRoot.py | python3 TreeCluster.py -m $mode -t $thresh > $clusters
for t in 8 9; do 
    python3 cluster_previous_time.py -c $clusters -tn $trans -t $t > "$clusters.time$t"
done
python3 cluster_growth_rates.py -c1 $time8 -c2 $time9 > $growth
python3 pick_by_cluster_growth.py -c $time9 -g $growth -n $num | cut -d'|' -f2 > $people
python3 cluster_efficacy.py -i $people -tn $trans -t 9

# clean up
rm -f $clusters $time8 $time9 $growth $people
