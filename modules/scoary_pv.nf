// Filter panaroo output (pvs) with scoary
process scoary_pv {
    publishDir  "${params.outdir}/2.genomic_features/pv_out/scoary_pv", mode: 'copy', overwrite: true
    cache 'lenient'

    input:
    path pv_coll_traitfile
    path pvs 

    output:
    path '*.results.csv' 

    script:
    """
    scoary -t $pv_coll_traitfile -g $pvs --no-time --collapse -p 1.0
    """
}
