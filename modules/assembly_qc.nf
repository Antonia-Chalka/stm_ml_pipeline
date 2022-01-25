process assembly_qc {
    cache 'lenient'

    input:
    path assembly

    output:
    path "${assembly.baseName}_qc/transposed_report.tsv"

    script:
    """
    quast.py  -o "${assembly.baseName}_qc" "${assembly}" 
    """
}
