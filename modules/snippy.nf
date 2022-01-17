process snippy {
    cache 'lenient'

    input:
    path assembly 
    path snp_ref_file

    output:
    path "./${assembly.baseName}/" 

    script:
    """
    snippy --prefix "${assembly.baseName}" --outdir "${assembly.baseName}" --ref  "${snp_ref_file}" --ctgs "${assembly}" --force
    """
}
