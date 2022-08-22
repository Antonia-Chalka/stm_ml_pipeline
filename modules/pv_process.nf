process pv_process {
    publishDir  "${params.outdir}/4.model/model_input", mode: 'copy', overwrite: true, pattern: "*.tsv"
    publishDir  "${params.outdir}/2.genomic_features/pv_out/scoary_pv", mode: 'copy', overwrite: true, pattern: "*.csv"

    cache 'lenient'

    input:
    path pv_process_script
    path pv_file_all
    path scoary_files
    path good_nonclonal_metadata
    
    output:
    path "pv_all.tsv", emit: pv_all
    path "*.tsv", emit: pv_inputs 

    script:
    """
    # Combine all scoary files into one & add headers
    for file in *.results.csv
    do
        sed '1d' \$file >>  scoary_all.csv
    done 
    sed -i "1s/^/\$(head -n1 \$file)\\n/" scoary_all.csv

    Rscript --vanilla $pv_process_script ${pv_file_all} scoary_all.csv ${good_nonclonal_metadata} $params.assembly_column $params.host_column $params.scoarycutoff
    """
}
