process snp_process_test {
    publishDir  "${params.outdir}/4.model/model_input", mode: 'copy', overwrite: true, pattern: "*.tsv"
    cache 'lenient'

    input:
    path snps_process_script
    path snp_core
    path snp_core_ref
    
    output:
    path "snp_abudance_all.tsv"

    script:
    """
    Rscript --vanilla $snps_process_script $snp_core $snp_core_ref
    """
}
