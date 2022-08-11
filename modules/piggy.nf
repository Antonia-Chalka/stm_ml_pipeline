// Run piggy (extracts intergenic regions)
process piggy {
    publishDir  "${params.outdir}/2.genomic_features/igr_out", mode: 'copy', overwrite: true
    cache 'lenient'

    input:
    path annotations // from annotation_piggy.collect()
    path roary_dir 

    output:
    path "piggy_out/IGR_presence_absence.csv", emit: piggy_csv
    path "piggy_out/IGR_presence_absence.Rtab", emit: piggy_rtab
    path "piggy_out/representative_clusters_merged.fasta", emit: piggy_ref
    path "piggy_out/cluster_intergenic_alignment_files/", emit: piggy_igr_aln_dir
    path "piggy_out/switched_region_alignment_files/", emit: piggy_switched_aln_dir
    path "piggy_out/cluster_IGR_divergences.csv", emit: cluster_IGR_divergences
    path "piggy_out/core_IGR_alignment.fasta", emit:core_IGR_alignment
    path "piggy_out/IGR_sequences.fasta", emit: IGR_sequences
    path "piggy_out/log.txt", emit: log
    path "piggy_out/representative_clusters_merged.fasta", emit: representative_clusters_merged
    path "piggy_out/roary_piggy_combined.tab", emit: roary_piggy_combined
    path "piggy_out/switched_region_divergences.csv", emit: switched_region_divergences

    script:
    """
    piggy -r "${roary_dir}" -i . -t $params.threads
    """
}
