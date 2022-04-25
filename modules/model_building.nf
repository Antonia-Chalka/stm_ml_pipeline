// Model generation - Host/Source Attribution
process model_building {
    publishDir "${params.outdir}/4.model/models", mode: 'copy', overwrite: true, pattern: '*.rds'
    publishDir "${params.outdir}/4.model/predictions", mode: 'copy', overwrite: true, pattern: '*.csv'
    cache 'lenient'

    input:
    path model_building_script
    path input

    output:
    path '*.rds', emit: model
    path '*.csv', emit: predictions

    script:
    """
    Rscript --vanilla $model_building_script $input
    """
}
