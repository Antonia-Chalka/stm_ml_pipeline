process clonal_filtering {
    publishDir  "${params.outdir}/clonal_detection", mode: 'copy', overwrite: true, pattern : "*.txt"
    cache 'lenient'
    
    input:
    path clusters
    path good_assemblies_list
    path metadata

    output:
    path "good_nonclonal_metadata.csv"

    script:
    """
    cut -f1 ${good_assemblies_list} > good_assemblies_all.lst
    cat *.list > cluster_all.txt

    # If no clonal clusters were detected, exit
    if ! [ -s "cluster_all.txt" ];then
        cat $metadata > good_nonclonal_metadata.csv 
        exit
    fi

    # Pick 1 random line from the input file (use as chsen seqs)
    for file in *.list; 
    do 
        head -\$((\${RANDOM} % `wc -l < "\$file"` + 1)) "\$file" | tail -1
    done > cluster_chosen.txt

    # take intersection of above filen
    grep -Fxvf cluster_chosen.txt cluster_all.txt > cluster_discarded.txt
    
    # generate list of all sequences without discarded
    grep -vf cluster_discarded.txt  good_assemblies_all.lst  > all_nonclonal.txt

    # Create a new metadata file of the filtered results & add headers
    awk -F',' 'NR==FNR{c[\$1]++;next};c[\$1] > 0' all_nonclonal.txt $metadata > good_nonclonal_metadata.csv 
    printf '%s\n' '0r !head -n 1 $metadata' x | ex good_nonclonal_metadata.csv 
    """
}
