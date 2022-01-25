// Process AMR Data for model generation
process amr_process {
    publishDir  "${params.outdir}/model_input", mode: 'copy', overwrite: true, pattern: '*.tsv'
    cache 'lenient'

    input:
    path amr_process_script
    path amr_file_all
    path good_nonclonal_metadata
    
    output:
    path "amr_gene_all.tsv", emit: amr_gene_all
    path "amr_gene_bps.tsv", emit: amr_gene_bps
    path "amr_gene_human.tsv", emit: amr_gene_human
    path "amr_class_all.tsv", emit: amr_class_all
    path "amr_class_bps.tsv", emit: amr_class_bps
    path "amr_class_human.tsv", emit: amr_class_human

    script:
    """
    Rscript --vanilla $amr_process_script ${amr_file_all} ${good_nonclonal_metadata} $params.assembly_column $params.host_column
    """
}
