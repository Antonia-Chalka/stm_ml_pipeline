process get_igr_fastas {
    publishDir  "${params.outdir}/4.model/model_ref", mode: 'copy', overwrite: true, pattern : "igr_filter.fasta"
    cache 'lenient'
    
    input:
    path igr_model_data
    path piggy_igr_ref_file

    output:
    path "igr_filter.fasta"
    path "igr_missing_headers.txt"

    script:
    """
    # Clean headers of representative igr file
    sed '/>/ s/_+.*//g' $piggy_igr_ref_file > cleaned.fasta

    # Get names of IGRs inputted in the model
    head -n 1 $igr_model_data | sed 's/\"//g' | tr -s ' \t' '\n' | head -n -1 > igr.list

    # Filter list of igrs
    seqtk subseq cleaned.fasta igr.list > igr_filter.fasta

    # Check if the filtered list contains any missing headers
    grep ">" igr_filter.fasta | cut -c 2- > filtered_headers.fasta
    comm -13 <(sort filtered_headers.fasta) <(sort igr.list) > igr_missing_headers.txt
    """
}
