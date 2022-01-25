//Concantenate all amr files
process amr_collect {
    publishDir  "${params.outdir}/amr_out", mode: 'copy', overwrite: true, pattern : "amr_all.tsv"
    cache 'lenient'

    input:
    path amr_files

    output:
    path "amr_all.tsv"

    script:
    """
    for file in *.tsv
    do
        name=\${file%_amr.tsv}
        sed -i '1d' \$file
        sed "s/\$/\t\$name/" "\$file"
    done > amr_all.tsv

    # Add headers
    sed -i.bak 1i"Name\tProtein identifier\tContig id\tStart\tStop\tStrand\tGene symbol\tSequence name\tScope\tElement.type\tElement subtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference sequence\t% Identity to reference sequence\tAlignment length\tAccession of closest sequence\tName of closest sequence\tHMM id\tHMM description\tFilename" amr_all.tsv
    """
}
