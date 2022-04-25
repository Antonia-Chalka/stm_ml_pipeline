// Process AMR Data for model generation
process amr_process {
    publishDir  "${params.outdir}/4.model/model_input", mode: 'copy', overwrite: true, pattern: '*.tsv'

    input:
    path amr_process_script
    path amr_file_all
    path good_nonclonal_metadata
    
    output:
    path '*.tsv', emit: amr_inputs

    script:
    """
    Rscript --vanilla $amr_process_script ${amr_file_all} ${good_nonclonal_metadata} $params.assembly_column $params.host_column
    """
}
