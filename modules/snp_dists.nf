// Get SNP distance of assemblies to use in nonclonal filtering
process snp_dists {
    cache 'lenient'

    input:
    path snp_core_aln 

    output:
    path 'snpdist_base.tsv' 

    script:
    """
    snp-dists -m $snp_core_aln > snpdist_base.tsv
    """
}
