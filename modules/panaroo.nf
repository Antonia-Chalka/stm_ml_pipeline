// Run panaroo 
process panaroo {
    publishDir  "${params.outdir}/panaroo_out", mode: 'copy', overwrite: true
    cache 'lenient'

    input:
    path annotations //from annotation_panaroo.collect()

    output:
    path "panaroo_out/gene_presence_absence_roary.csv", emit: pv_csv
    path "panaroo_out/gene_presence_absence.Rtab", emit: pv_rtab
    path "panaroo_out/pan_genome_reference.fa", emit: pv_ref
    path "roary_out/", emit: roary_dir

    script:
    """
    panaroo -i ${annotations} -o panaroo_out/ --clean-mode moderate -t $params.threads --remove-invalid-genes --merge_paralogs
    mkdir ./roary_out
    cp  ./panaroo_out/gene_presence_absence_roary.csv ./roary_out/gene_presence_absence.csv
    """
}
