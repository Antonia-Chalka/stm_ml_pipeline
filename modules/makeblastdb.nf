process makeblastdb {
    publishDir  "${params.outdir}/4.model/model_ref", mode: 'copy', overwrite: true, pattern : "./blastdb/"
    cache 'lenient'
    
    input:
    path assemblies

    output:
    path "blastdb/", emit: blastdb_dir

    script:
    """
    # Append filename at start of fasta headers & concantenate
    for f in *.${params.fileextension}; 
    do 
        sed "s|^>|>\${f} |g" "\${f}"; 
    done > all_assemblies.fa

    mkdir .blastdb/

    makeblastdb -in all_assemblies.fa -dbtype nucl  -out ./blastdb/assemblyblastdb
    
    rm -f all_assemblies.fasta
    """
}
