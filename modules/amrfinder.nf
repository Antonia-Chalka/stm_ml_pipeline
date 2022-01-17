// Run Amrfinder using all possible inputs (assembly, ggff, and trans portein files)
process amrfinder {
    cache 'lenient'

    input:
    path annotation 
    path trans_protein 
    path nucleotide 

    output:
    path "${annotation.baseName}_amr.tsv" 

    script:
    """
    # Convert the prokka generated gff into a format amrfinder accepts
    perl -pe '/^##FASTA/ && exit; s/(\\W)Name=/\$1OldName=/i; s/ID=([^;]+)/ID=\$1;Name=\$1/' $annotation  > amrfinder.gff

    amrfinder -p ${trans_protein} -g amrfinder.gff -n ${nucleotide} -O ${params.amr_species} -o "${annotation.baseName}_amr.tsv" --name --plus
    """
}
