process get_pv_fastas {
    publishDir  "${params.outdir}/4.model/model_ref", mode: 'copy', overwrite: true, pattern : "pv_filter.fasta"
    cache 'lenient'
    
    input:
    path pv_model_data
    path pv_clusters

    output:
    path "pv_filter.fasta"

    script:
    """
    # Get names of pvs inputted in the model
    head -n 1 $pv_model_data | sed 's/\"//g' | tr -s ' \t' '\n' | head -n -1 > pv.list

    # Filter list of pvs
    seqtk subseq $pv_clusters pv.list > pv_filter.fasta

    # Check if the filtered list contains any missing headers
    grep ">" pv_filter.fasta | cut -c 2- > filtered_headers.fasta
    comm -13 <(sort filtered_headers.fasta) <(sort pv.list) > missing_headers.txt
    """
}
