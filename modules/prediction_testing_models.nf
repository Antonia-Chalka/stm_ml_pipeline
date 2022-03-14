// Model generation - Host/Source Attribution
process prediction_testing_models {
    publishDir "${params.outdir}/4.model/predictions", mode: 'copy', overwrite: true, pattern: '*.tsv'

    input:
    path amr_class_all
    path amr_gene_all
    path pv_blast
    path igr_blast
    path snp_abudance_all
    path model_dir
    path model_testing_script

    output:
    path '*.tsv', emit: predictions

    script:
    """
    Rscript --vanilla $model_testing_script $amr_class_all $amr_gene_all $pv_blast $igr_blast $snp_abudance_all
    """
}
