process igr_process {
    publishDir  "${params.outdir}/4.model/model_input", mode: 'copy', overwrite: true, pattern: "*.tsv"
    cache 'lenient'

    input:
    path igr_process_script
    path igr_file_all
    path scoary_files
    path good_nonclonal_metadata
    
    output:
    path "igr_all.tsv", emit: igr_all
    path "*.tsv", emit: igr_inputs

    script:
    """
    # Combine all scoary files into one & add headers
    for file in *.results.csv
    do
        sed '1d' \$file >>  scoary_all.csv
    done 
    sed -i "1s/^/\$(head -n1 \$file)\\n/" scoary_all.csv

    Rscript --vanilla $igr_process_script ${igr_file_all} scoary_all.csv ${good_nonclonal_metadata} $params.assembly_column $params.host_column $params.scoarycutoff
    """
}
