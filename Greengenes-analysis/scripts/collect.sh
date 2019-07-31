for alg in max sum single-linkage;  do for c in `ls ${alg}`; do printf "patristic\t%s\tTREECLUSTER-${alg}\t" $c; cat $alg/$c/weighted_mean.txt| tr '\n' '\t'; /home/balaban/clusterproj/scripts/uc_count.sh $alg/$c/clusters_formatted.txt | awk '{print $3}' ; printf "hamming\t%s\tTREECLUSTER-${alg}\t" $c; cat ${alg}/$c/weighted_mean_hamming.txt | tr '\n' '\t'; /home/balaban/clusterproj/scripts/uc_count.sh ${alg}/$c/clusters_formatted.txt | awk '{print $3}'; done; done 

for alg in gg; do
    for c in `ls ${alg}`; do
        printf "patristic\t%s\tGG(UCLUST)\t" $c; cat $alg/${c}/weighted_mean.txt | tr '\n' '\t'; /home/balaban/clusterproj/scripts/uc_count.sh $alg/${c}/clusters_formatted.txt | awk '{print $3}' ; printf "hamming\t%s\tGG(UCLUST)\t" ${c}; cat $alg/${c}/weighted_mean_hamming.txt | tr '\n' '\t'; /home/balaban/clusterproj/scripts/uc_count.sh $alg/${c}/clusters_formatted.txt | awk '{print $3}'
    done
 done
