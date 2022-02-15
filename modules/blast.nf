process blast {
    publishDir  "${params.outdir}/blast_results", mode: 'copy', overwrite: true, pattern : "*_results.tsv"
    cache 'lenient'
    
    input:
    path blastdb_dir
    path query_seq

    output:
    path "${query_seq}_results.tsv", emit: blast_results

    script:
    """
    blastn -query ${query_seq} -db ./${blastdb_dir}/assemblyblastdb -outfmt 6 -num_threads ${params.blast_threads} -out "${query_seq}_results.tsv"
    """
}
