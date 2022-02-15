// Filter assemblies based on quast-derived metrics
process printqc {
    publishDir  "${params.outdir}/1.assembly_quality/good_assemblies/", mode: 'copy', overwrite: true, pattern : "*.${fileextension}"
    publishDir  "${params.outdir}/1.assembly_quality/", mode: 'copy', overwrite: true, pattern : 'good_noext_metadata.csv'

    cache 'lenient'

    input:
    path qcs
    path metadata_file

    output:
    path "*.${params.fileextension}", emit: good_assemblies
    path 'good_noext_metadata.csv', emit: good_metadata
    path 'pass_qc.tsv', emit: good_assemblies_list

    script:
    """
    # Create list of assemblies that pass qc & Copy them
    awk -F "\t" '{ if((\$2 < $params.ctg_count) && (\$8 > $params.as_ln_lwr && \$8 < $params.as_ln_upr) && (\$15 > $params.largest_ctg) && (\$17 > $params.gc_lwr && \$17 < $params.gc_upr) && (\$18 > $params.n50)) { print } }' $qcs > pass_qc.tsv
    for assembly_name in `cut -f1 pass_qc.tsv`
    do 
        ln -s ${params.assemblypath}/\${assembly_name}.${params.fileextension} \${assembly_name}.${params.fileextension}
    done

    # Create list from good assemblies & create a new metadata file of the filtered results
    ls *.${params.fileextension} > good_assemblies.txt
    awk -F',' 'NR==FNR{c[\$1]++;next};c[\$1] > 0' good_assemblies.txt $metadata_file > good_metadata.csv
    # Add headers to new metadata file by copying the first line of the orginal metadata file
    printf '%s\n' '0r !head -n 1 $metadata_file' x | ex good_metadata.csv 

    # Remove extensions from metadata file
    awk '{ gsub(/\\.${params.fileextension}/,"", \$1); print }' good_metadata.csv  > good_noext_metadata.csv
    """
}
