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
params.prokka_ref = "$projectDir/data/stm_proteinref.fasta"
prokka_ref_file = file(params.prokka_ref)

// AMRfinder parameters
params.amr_species='Salmonella'

// Panaroo parameters
params.panaroo_mode = 'moderate'

// Snippy parameters
params.snp_ref = "$projectDir/data/stm_sl1344.fasta"
snp_ref_file = file(params.snp_ref)

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
        ln -s ${params.assemblypath}/\${assembly_name}.fasta ./\${assembly_name}.fasta 
    done
    """
}

good_assemblies.into {assemblies_prokka; assemblies_snippy}

// Annotate via prokka. Reference protein file required
process prokka_annotation {
    publishDir  "${params.outdir}/annotation", mode: 'copy', overwrite: true

    input:
    file good_qc_assembly from assemblies_prokka.flatten()
    file prokka_ref_file

    output:
    file "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.gff" into annotation_ch
    file "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.faa" into trans_protein_ch
    file "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.fna" into nucleotide_ch

    script:
    """
    prokka --prefix "${good_qc_assembly.baseName}" --force --proteins ${prokka_ref_file} --centre X --compliant ${good_qc_assembly}
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

process piggy {
    publishDir  "${params.outdir}", mode: 'copy', overwrite: true

    input:
    file annotations from annotation_piggy.collect()
    file roary_dir from roary_dir_ch

    output:
    file "./piggy_out/IGR_presence_absence.csv" into piggy_igr_ch

    script:
    """
    piggy -r "${roary_dir}" -i .
    """
}

// scoary filtering - need to make host score data (Assembly,bovine,human,poultry,swine) via r script

// panaroo - do qc script as well?
process panaroo {
    input:
    file annotations from annotation_panaroo.collect()

    output:
    file "panaroo_out/" into panaroo_dir_ch

    script:
    """
    panaroo -i ${annotations} -o panaroo_out/ --clean-mode moderate
    """
}

// scoary

// snippy
process snippy {
    input:
    file assembly from assemblies_snippy.flatten()
    file snp_ref_file

    output:
    file "./${assembly.baseName}/" into snippy_folders_ch

    script:
    """
    snippy --prefix "${assembly.baseName}" --outdir "${assembly.baseName}" --ref  "${snp_ref_file}" --ctgs "${assembly}" --force
    """
}

process snippy_core {
    input:
    file snippy_folders from snippy_folders_ch.collect()
    file snp_ref_file

    output:
    file 'core.tab' into snps_ch

    script:
    """
    snippy-core --prefix core --ref "${snp_ref_file}" ${snippy_folders}
    """
}