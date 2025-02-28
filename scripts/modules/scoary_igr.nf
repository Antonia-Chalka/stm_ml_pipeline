// Filter piggy output (igrs) with scoary
process scoary_igr {
    publishDir  "${params.outdir}/2.genomic_features/igr_out/scoary_igr", mode: 'copy', overwrite: true
    cache 'lenient'

    input:
    path igr_coll_traitfile
    path igrs 

    output:
    path '*.results.csv' 

    script:
    """
    scoary -t $igr_coll_traitfile -g $igrs --no-time --collapse -p 1.0
    """
}
