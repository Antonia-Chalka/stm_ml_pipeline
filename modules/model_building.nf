// Model generation - Host/Source Attribution
process model_building {
    publishDir "${params.outdir}/4.model/models", mode: 'copy', overwrite: true, pattern: '*.rds'
    publishDir "${params.outdir}/4.model/predictions", mode: 'copy', overwrite: true, pattern: '*.csv'
    cache 'lenient'

    input:
    path amr_class_all
    path amr_class_bps
    path amr_gene_all
    path amr_gene_bps
    path pv_all
    path pv_bps
    path igr_all
    path igr_bps
    path snp_abudance_all
    path snp_abudance_bps
    path model_building_script

    output:
    path '*.rds', emit: models
    path '*.csv', emit: predictions

    script:
    """
    Rscript --vanilla $model_building_script $amr_class_all $amr_class_bps $amr_gene_all $amr_gene_bps $pv_all $pv_bps $igr_all $igr_bps $snp_abudance_all $snp_abudance_bps
    """
}
