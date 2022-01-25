process snp_process {
    publishDir  "${params.outdir}/model_input", mode: 'copy', overwrite: true, pattern: "*.tsv"
    cache 'lenient'

    input:
    path snps_process_script
    path snp_core
    path good_nonclonal_metadata
    
    output:
    path "snp_abudance_all.tsv", emit: snp_abudance_all
    path "snp_abudance_bps.tsv", emit: snp_abudance_bps
    path "snp_abudance_human.tsv", emit: snp_abudance_human

    script:
    """
    Rscript --vanilla $snps_process_script $snp_core $good_nonclonal_metadata $params.assembly_column $params.host_column
    """
}
