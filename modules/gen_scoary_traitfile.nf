// scoary filtering - make traitfile(Assembly,[Hosts]]) via r script
process gen_scoary_traitfile {
    publishDir  "${params.outdir}/2.genomic_features/", mode: 'copy', overwrite: true
    cache 'lenient'

    input:
    path metadata
    path scoary_datagen_file

    output:
    path "scoary_traitfile.csv" 

    script:
    if (scoary_datagen_file=="$projectDir/data/scoary_generate_tabfile.R")
    """
    Rscript --vanilla $scoary_datagen_file $metadata $params.assembly_column $params.host_column
    """
    else
    """
    Rscript --vanilla $scoary_datagen_file $metadata $params.assembly_column $params.host_column $params.scoary_threshold
    """

}
