// Run piggy (extracts intergenic regions)
process piggy {
    publishDir  "${params.outdir}", mode: 'copy', overwrite: true
    cache 'lenient'

    input:
    path annotations // from annotation_piggy.collect()
    path roary_dir 

    output:
    path "./piggy_out/IGR_presence_absence.csv", emit: piggy_csv
    path "./piggy_out/IGR_presence_absence.Rtab", emit: piggy_rtab
    path "./piggy_out/representative_clusters_merged.fasta", emit: piggy_ref

    script:
    """
    piggy -r "${roary_dir}" -i . -t $params.threads
    """
}
