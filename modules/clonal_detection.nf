// Filter for snps
// HACK If no clusters are detected, no list file is produced. so a blank one has been used for now
process clonal_detection {
    publishDir  "${params.outdir}/clonal_detection", mode: 'copy', overwrite: true, pattern : "*.png"
    cache 'lenient'

    input:
    path clonal_detection_script
    path metadata
    path snp_dist_file 

    output:
    path "*.list"

    script:
    """
    touch empty.list
    Rscript --vanilla $clonal_detection_script $snp_dist_file $metadata $params.assembly_column $params.host_column $params.year_collected $params.region_column $params.snp_dist_threshold
    """
}