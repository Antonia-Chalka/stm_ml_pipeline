// Run panaroo 
process panaroo {
    publishDir  "${params.outdir}/2.genomic_features/pv_out", mode: 'copy', overwrite: true
    cache 'lenient'

    input:
    path annotations //from annotation_panaroo.collect()

    output:
    path "panaroo_out/gene_presence_absence_roary.csv", emit: pv_csv
    path "panaroo_out/gene_presence_absence.Rtab", emit: pv_rtab
    path "panaroo_out/pan_genome_reference.fa", emit: pv_ref
    path "roary_out/", emit: roary_dir
    path "panaroo_out/combined_DNA_CDS.fasta", emit: combined_DNA_CDS
    path "panaroo_out/combined_protein_cdhit_out.txt", emit: combined_protein_cdhit_out
    path "panaroo_out/combined_protein_cdhit_out.txt.clstr", emit: combined_protein_cdhit_out_clstr
    path "panaroo_out/combined_protein_CDS.fasta", emit: combined_protein_CDS
    path "panaroo_out/final_graph.gml", emit: final_graph
    path "panaroo_out/gene_data.csv", emit: gene_data
    path "panaroo_out/gene_presence_absence.csv", emit: gene_presence_absence
    path "panaroo_out/pre_filt_graph.gml", emit: pre_filt_graph
    path "panaroo_out/struct_presence_absence.Rtab", emit: struct_presence_absence
    path "panaroo_out/summary_statistics.txt", emit: summary_statistics

    script:
    """
    mkdir ./roary_out
    
    panaroo -i ${annotations} -o panaroo_out/ --clean-mode moderate -t $params.threads --remove-invalid-genes --merge_paralogs
    
    sed -e 's/,/","/g' ./panaroo_out/gene_presence_absence_roary.csv > ./roary_out/gene_presence_absence.csv
    """
}
