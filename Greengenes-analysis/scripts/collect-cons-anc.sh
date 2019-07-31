for alg in max;  do 
    for c in `ls ${alg}`; do 
        k=`/home/balaban/clusterproj/scripts/uc_count.sh $alg/$c/clusters_formatted.txt | awk '{print $3}'`
        cat $alg/$c/consensus_all.txt | sed 's/_/\t/' | sed "s/^/$c\tTC-max\t/" | sed "s/$/\t$k/" | sed "s/TC-max\tconsensus/TC-max-consensus/"
        cat $alg/$c/ancestral_all.txt | sed 's/_/\t/' | sed "s/^/$c\tTC-max\t/" | sed "s/$/\t$k/" | sed "s/TC-max\tancestral/TC-max-ancestral/"
    done; 
done 

for alg in gg;  do 
    for c in `ls ${alg}`; do 
        k=`/home/balaban/clusterproj/scripts/uc_count.sh $alg/$c/clusters_formatted.txt | awk '{print $3}'`
        cat $alg/$c/consensus_all.txt | sed 's/_/\t/' | sed "s/^/$c\tGG\t/" | sed "s/$/\t$k/" | sed "s/GG\tconsensus/GG-centroid/"
    done; 
done 


#for alg in gg; do
#    for c in `ls ${alg}`; do
#        printf "patristic\t%s\tGG(UCLUST)\t" $c; cat $alg/${c}/weighted_mean.txt | tr '\n' '\t'; /home/balaban/clusterproj/scripts/uc_count.sh $alg/${c}/clusters_formatted.txt | awk '{print $3}' ; printf "hamming\t%s\tGG(UCLUST)\t" ${c}; cat $alg/${c}/weighted_mean_hamming.txt | tr '\n' '\t'; /home/balaban/clusterproj/scripts/uc_count.sh $alg/${c}/clusters_formatted.txt | awk '{print $3}'
#    done
# done
