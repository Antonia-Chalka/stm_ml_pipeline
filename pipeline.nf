#!/usr/bin/env nextflow

// Required params - defaults to testing data
// HACK THE DIRECTORIES HAVE TO BE ABSOLUTE ONES
params.outdir = "$projectDir/out"
params.assemblypath = "$projectDir/input_test/"
params.hostdata = "$projectDir/input_test/metadata.csv"

params.assembly_column="Filename"
params.host_column="Source.Host"
params.fileextension="fasta"

metadata_file = file(params.hostdata)
assemblies = Channel.fromPath("${params.assemblypath}/*.${params.fileextension}")

// Assembly quality thresholds
params.as_ln_upr = 6000000
params.as_ln_lwr = 4000000
params.ctg_count = 500
params.largest_ctg = 100000
params.n50 = 50000
params.gc_upr = 54
params.gc_lwr = 50

// Computing parameters
params.threads = 5

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

// Scoary hostfile script
scoary_datagen_file=file("$projectDir/data/scoary_generate_tabfile.R")

// R scripts for model input processing
amr_process_file=file("$projectDir/data/input_amr.R")
pv_process_file=file("$projectDir/data/input_pv.R")
igr_process_file=file("$projectDir/data/input_igr.R")
snps_process_file=file("$projectDir/data/input_snps.R")

// Run Quast
process assembly_qc {
    publishDir  "${params.outdir}/qc_report", mode: 'copy', overwrite: true

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
    publishDir  "${params.outdir}/good_assemblies", mode: 'copy', overwrite: true, pattern : "${params.fileextension}"
    input:
    file qcs from quast_ch.collectFile(keepHeader:true, skip:1)
    file metadata_file

    output:
    file "*.${params.fileextension}" into good_assemblies
    file 'good_metadata.csv' into good_metadata

    script:
    """
    # Create list of assemblies that pass qc & Copy them
    awk -F "\t" '{ if((\$2 < $params.ctg_count) && (\$8 > $params.as_ln_lwr && \$8 < $params.as_ln_upr) && (\$15 > $params.largest_ctg) && (\$17 > $params.gc_lwr && \$17 < $params.gc_upr) && (\$18 > $params.n50)) { print } }' $qcs > pass_qc.tsv
    for assembly_name in `cut -f1 pass_qc.tsv`
    do 
        ln -s ${params.assemblypath}/\${assembly_name} \${assembly_name} 
    done

    # Create list from good assemblies & create a new metadata file of the filtered results
    ls *.${params.fileextension} > good_assemblies.txt
    awk -F',' 'NR==FNR{c[\$1]++;next};c[\$1] > 0' good_assemblies.txt $metadata_file > good_metadata.csv
    # Add headers to new metadata file by copying the first line of the orginal metadata file
    printf '%s\n' '0r !head -n 1 $metadata_file' x | ex good_metadata.csv 
    """
}

good_assemblies.into {assemblies_prokka; assemblies_snippy}
good_metadata.into {scoary_metadata; amr_process_metadata; model_metadata}

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


// roary
process roary {

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
    piggy -r "${roary_dir}" -i . -t $params.threads
    """
}

// scoary filtering - make traitfile(Assembly,[Hosts]]) via r script
// TODO move scoary datagen script elsewehere?
process gen_scoary_traitfile {
    publishDir  "${params.outdir}/", mode: 'copy', overwrite: true

    input:
    file scoary_metadata
    file scoary_datagen_file

    output:
    file "./scoary_traitfile.csv" into scoary_traitfile_ch

    script:
    """
    Rscript --vanilla $scoary_datagen_file $scoary_metadata $params.assembly_column $params.host_column $params.fileextension
    """
}

scoary_traitfile_ch.into { pv_coll_traitfile; igr_coll_traitfile }

// scoary pv TODO
process scoary_igr {
    publishDir  "${params.outdir}/scoary_igr", mode: 'copy', overwrite: true

    input:
    file igr_coll_traitfile
    file igrs from piggy_igr_ch

    output:
    file '*.results.csv' into igr_scores

    script:
    """
    scoary -t $igr_coll_traitfile -g $igrs --no-time --collapse -p 1.0
    """
}

// panaroo - TODO do qc script as well?
process panaroo {
    publishDir  "${params.outdir}/panaroo_out", mode: 'copy', overwrite: true

    input:
    file annotations from annotation_panaroo.collect()

    output:
    file "panaroo_out/gene_presence_absence_roary.csv" into panaroo_pres_abs_ch

    script:
    """
    panaroo -i ${annotations} -o panaroo_out/ --clean-mode moderate -t $params.threads
    """
}

// scoary 
process scoary_pv {
    publishDir  "${params.outdir}/scoary_pv", mode: 'copy', overwrite: true

    input:
    file pv_coll_traitfile
    file pvs from panaroo_pres_abs_ch

    output:
    file '*.results.csv' into pv_scores

    script:
    """
    scoary -t $pv_coll_traitfile -g $pvs --no-time --collapse -p 1.0
    """
}

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
    publishDir  "${params.outdir}/snippy_core_out", mode: 'copy', overwrite: true

    input:
    file snippy_folders from snippy_folders_ch.collect()
    file snp_ref_file

    output:
    file 'core.tab' into snps_ch
    file 'core.aln' into snp_core_aln_ch

    script:
    """
    snippy-core --prefix core --ref "${snp_ref_file}" ${snippy_folders}
    """
}

// TODO Add Data prep for model generation
/*
process amr_process {
    publishDir  "${params.outdir}/model_input", mode: 'copy', overwrite: true

    input:
    file amr_process_metadata
    file amr_file_all from amr_single_ch.collectFile(keepHeader:true, skip:1, storeDir:"${params.outdir}/amr_out")

    output:
    file "*.tsv" into amr_model_inputs_ch

    script:
    """
    cat $amr_file_all
    Rscript --vanilla ${amr_file_all} ${amr_process_metadata} $params.assembly_column $params.host_column
    """
}
*/

// .collectFile(keepHeader:true, skip:1, name:"amr_all.tsv", storeDir:"${params.outdir}/amr_out")
// TODO Add basic model generation
// # Rscript --vanilla ${amr_file_all} $params.assembly_column $params.host_column

// TODO ADD SNP-DIST for Nonclonal filtering - but, need to have Region & Year -NOT URRENTLY IMPLEMENTED
process snp_dists {

    input:
    file snp_core_aln from snp_core_aln_ch

    output:
    file 'snpdist_base.tsv' into snp_dist_ch

    script:
    """
    snp-dists -m $snp_core_aln > snpdist_base.tsv
    """
}