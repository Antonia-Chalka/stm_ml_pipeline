process snippy_core {
    publishDir  "${params.outdir}/2.genomic_features/snp_core_out", mode: 'copy', overwrite: true
    cache 'lenient'

    input:
    path snippy_folders
    path snp_ref_file

    output:
    path 'core.tab', emit: core_snps
    path 'core.aln', emit: aligned_snps
    path 'core.full.aln', emit: core_full_snps
    path 'clean.full.aln', emit: core_full_clean_snps
    path 'core.ref.tab', emit: core_ref_tab

    script:
    """
    snippy-core --prefix core --ref "${snp_ref_file}" ${snippy_folders}
    
    # produce 'cleaned' full aln file to use in phylogeny
    snippy-clean_full_aln core.full.aln > clean.full.aln
    
    # Get the SNP bases for the reference to use when testing new sequences
    cut -f 1-3 core.tab > core.ref.tab
    """
}
