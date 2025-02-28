// Process AMR Data for model generation
process amr_process_test {
    publishDir  "${params.outdir}/4.model/model_input", mode: 'copy', overwrite: true, pattern: '*.tsv'
    cache 'lenient'

    input:
    path amr_process_script
    path amr_file_all
    
    output:
    path "amr_gene_all.tsv", emit: amr_gene_all
    path "amr_class_all.tsv", emit: amr_class_all

    script:
    """
    Rscript --vanilla $amr_process_script ${amr_file_all}
    """
}
