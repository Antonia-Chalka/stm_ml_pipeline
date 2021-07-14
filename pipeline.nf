#!/usr/bin/env nextflow

// Required params - defaults to testing data
params.hostdata = "$projectDir/input_test/metadata.csv"
params.outdir = "$projectDir/out"
params.assemblypath = "$projectDir/input_test/"
assemblies = Channel.fromPath("${params.assemblypath}/*.fasta")

// Assembly quality thresholds
params.as_ln_upr = 6000000
params.as_ln_lwr = 4000000
params.ctg_count = 500
params.largest_ctg = 100000
params.n50 = 50000
params.gc_upr = 54
params.gc_lwr = 50

// Prokka reference file - required
params.prokka_ref="$projectDir/data/stm_proteinref.fasta"

// AMRfinder parameters
params.amr_species='Salmonella'

// Run Quast
process assembly_qc {
    
    input:
    file assembly from assemblies

    output:
    file "${assembly.baseName}_qc/transposed_report.tsv" into quast_ch

    script:
    """
    quast.py  -o "${assembly.baseName}_qc" "${assembly}" 
    """
}

// Filter assemblies based on quast-derived metrics
process printqc {
    publishDir  "${params.outdir}/qc_report", mode: 'copy', overwrite: true, pattern : 'pass_qc.tsv'
    publishDir  "${params.outdir}/good_assemblies", mode: 'copy', overwrite: true, pattern : '*.fasta'

    input:
    file qcs from quast_ch.collectFile(keepHeader:true, skip:1)

    output:
    file '*.fasta' into good_assemblies

    script:
    """
    awk -F "\t" '{ if((\$2 < $params.ctg_count) && (\$8 > $params.as_ln_lwr && \$8 < $params.as_ln_upr) && (\$15 > $params.largest_ctg) && (\$17 > $params.gc_lwr && \$17 < $params.gc_upr) && (\$18 > $params.n50)) { print } }' $qcs > pass_qc.tsv

    for assembly_name in `cut -f1 pass_qc.tsv`
    do 
        echo \$assembly_name
        ln -s ${params.assemblypath}/\${assembly_name}.fasta ./\${assembly_name}.fasta 
    done
    """
}

// Annotate via prokka. Reference protein file required
process prokka_annotation {
    publishDir  "${params.outdir}/annotation", mode: 'copy', overwrite: true

    input:
    file good_qc_assembly from good_assemblies.flatten()

    output:
    file "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.gff" into annotation_ch
    file "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.faa" into trans_protein_ch
    file "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.fna" into nucleotide_ch

    script:
    """
    prokka --prefix "${good_qc_assembly.baseName}" --force --proteins ${params.prokka_ref} --centre X --compliant ${good_qc_assembly}
    """
}

annotation_ch.into { annotation_amrfinder; annotation_roary; annotation_piggy; annotation_panaroo }

process amrfinder {

    input:
    file annotation from annotation_amrfinder
    file trans_protein from trans_protein_ch
    file nucleotide from nucleotide_ch

    output:
    file "${annotation.baseName}_amr.tsv" into amr_single_ch

    script:
    """
    perl -pe '/^##FASTA/ && exit; s/(\\W)Name=/\$1OldName=/i; s/ID=([^;]+)/ID=\$1;Name=\$1/' $annotation  > amrfinder.gff

    amrfinder -p ${trans_protein} -g amrfinder.gff -n ${nucleotide} -O ${params.amr_species} -o "${annotation.baseName}_amr.tsv" --name --plus
    """
}

// collate amr files
amr_single_ch
    .collectFile(keepHeader:true, skip:1, name:"amr_all.tsv", storeDir:"${params.outdir}/amr_out")

// roary
process roary {
    publishDir  "${params.outdir}", mode: 'copy', overwrite: true

    input:
    file annotations from annotation_roary.collect()

    output:
    file "roary_out/" into roary_dir_ch

    script:
    """
    roary -v -s -f "roary_out" ${annotations}
    """
}
