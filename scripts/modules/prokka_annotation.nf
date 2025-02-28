// Annotate via prokka. Reference protein file required
process prokka_annotation {
    publishDir  "${params.outdir}/2.genomic_features/annotations", mode: 'copy', overwrite: true
    cache 'lenient'

    input:
    path good_qc_assembly 
    path prokka_ref_file

    output:
    path "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.gff", emit:annotation
    path "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.faa", emit:transl_protein
    path "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.fna", emit:nucleotide
    path "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.err", emit: err
    path "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.ffn", emit: ffn
    path "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.fsa", emit: fsa
    path "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.gbk", emit: gbk
    path "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.log", emit: log
    path "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.sqn", emit: sqn
    path "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.tbl", emit: tbl
    path "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.tsv", emit: tsv
    path "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.txt", emit: txt

    script:
    """
    prokka --prefix "${good_qc_assembly.baseName}" --force --proteins ${prokka_ref_file} --centre X --compliant ${good_qc_assembly} 
    """
}
