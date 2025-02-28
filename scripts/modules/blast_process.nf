process blast_process {
    publishDir  "${params.outdir}/4.model/model_input", mode: 'copy', overwrite: true, pattern: "*_blast.tsv"
    cache 'lenient'

    input:
    path blast_process_script
    path blast_input
    val outfile_prefix
    
    output:
    path "*_blast.tsv"

    script:
    """
    Rscript --vanilla $blast_process_script $blast_input $outfile_prefix
    """
}
