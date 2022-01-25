process snippy_core {
    publishDir  "${params.outdir}/snippy_core_out", mode: 'copy', overwrite: true
    cache 'lenient'

    input:
    path snippy_folders
    path snp_ref_file

    output:
    path 'core.tab', emit: core_snps
    path 'core.aln', emit: aligned_snps

    script:
    """
    snippy-core --prefix core --ref "${snp_ref_file}" ${snippy_folders}
    """
}
