process snp_process {
    publishDir  "${params.outdir}/4.model/model_input", mode: 'copy', overwrite: true, pattern: "*.tsv"
    cache 'lenient'

    input:
    path snps_process_script
    path snp_core
    path good_nonclonal_metadata
    
    output:
    path "*.tsv", emit: snp_inputs

    script:
    """
    Rscript --vanilla $snps_process_script $snp_core $good_nonclonal_metadata $params.assembly_column $params.host_column
    """
}
