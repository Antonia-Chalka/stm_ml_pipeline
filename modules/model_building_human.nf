// Model generation - Human scoring models
process model_building_human {
    publishDir "${params.outdir}/4.model/models", mode: 'copy', overwrite: true, pattern: '*.rds'
    publishDir "${params.outdir}/4.model/predictions", mode: 'copy', overwrite: true, pattern: '*.csv'
    publishDir "${params.outdir}/4.model/plots", mode: 'copy', overwrite: true, pattern: '*.png'

    input:
    path amr_class_human
    path amr_gene_human
    path pv_human
    path igr_human
    path snp_abudance_human
    path model_building_human_script

    output:
    path '*.rds', emit: models
    path '*.csv', emit: predictions
    path '*.png', emit: plots

    script:
    """
    Rscript --vanilla $model_building_human_script $amr_class_human $amr_gene_human $pv_human $igr_human $snp_abudance_human
    """
}
