#!/usr/bin/env nextflow
nextflow.enable.dsl=2

///////////////////////// Parameters /////////////////////////////////////////////////////////////////
// HACK Seeting it up as anything non-fasta messes up the quast output
fileextension="fasta"

// Use inputted prokka & snippy reference files
prokka_ref_file = file(params.prokka_ref)
snp_ref_file = file(params.snp_ref)

// R script to generate scoary traitfile
scoary_datagen_file=file("$projectDir/data/scoary_generate_tabfile.R")

// R script for clonal detection
clonal_detection_script=file("$projectDir/data/filter_script.R")

// R scripts for model input processing
amr_process_script=file("$projectDir/data/input_amr.R")
pv_process_script=file("$projectDir/data/input_pv.R")
igr_process_script=file("$projectDir/data/input_igr.R")
snps_process_script=file("$projectDir/data/input_snps.R")

// R scripts for model building
model_building_script=file("$projectDir/data/model_building.R")
model_building_human_script=file("$projectDir/data/model_building_human.R")


////////////////////// HELP Section //////////////////////////////////////////////
// TODO Change help message 
def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:

        --option                       Description/Notes [default value]

    Mandatory arguments:
         --assemblypath                Directory of your fasta files (full path required) [input_test]
         --hostdata                    CSV file containing assembly filename, host, year and region [input_test/metadata.csv]

    Optional arguments:
        --outdir                       Output directory for models & other data [./out]
        --snp_dist_threshold           SNP difference used to detect clonal clusters. Used in conjuction with region & collection year metadata. [10]
    
    Hostdata file options:
        --assembly_column              Column name of your assemblies. Must contain extension eg myassembly.fasta ["Filename"]
        --host_column                  Column name of your hosts ["Source.Host"]
        --region_column                Column of region of origin. Used to detect clonal clusters. Can be empty throughout whole dataset. ["Region"]
        --year_collection              Column of year each assembly was collected. Used to detect clonal clusters. Can be empty. ["Year"]

    Assembly Quality Thresholds:
        --as_ln_upr                    Maximum accepted assembly length [6000000]
        --as_ln_lwr                    Minumum accepted assembly length [4000000]
        --ctg_count                    Minimum accepted number of contigs [500]
        --largest_ctg                  Minimum accepted length of largest contig [100000] 
        --n50                          Minimum accepted n50 [50000]
        --gc_upr                       Minimum accepted GC % [54]
        --gc_lwr                       Maximum accepted GC % [50]

    Tool-specific arguments:
        --prokka_ref                   'Trusted' protein file for prokka (prokka --proteins) [./data/stm_proteinref.fasta]
        --amr_species                  Assembly species for amrfinder (amrfinder -O) ["Salmonella"]
        --panaroo_mode                 Panaroo assembly filtering mode (panaroo --clean-mode) ["moderate"]
        --snp_ref                      Reference file for snippy (snippy --ref) [/data/stm_sl1344.fasta]
        --threads                      Num of threads to use for panaroo & piggy (anaroo -t & piggy -t) [10]
        """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

////////////////////// Modules ////////////////////////////////////////////////////////////////////////
// Run Quast
process assembly_qc {
    publishDir  "${params.outdir}/qc_report", mode: 'copy', overwrite: true

    input:
    path assembly

    output:
    path "${assembly.baseName}_qc/transposed_report.tsv"

    script:
    """
    quast.py  -o "${assembly.baseName}_qc" "${assembly}" 
    """
}

// Filter assemblies based on quast-derived metrics
process printqc {
    publishDir  "${params.outdir}/good_assemblies", mode: 'copy', overwrite: true, pattern : "*.${fileextension}"
    input:
    path qcs
    path metadata_file

    output:
    path "*.${fileextension}", emit: good_assemblies
    path 'good_noext_metadata.csv', emit: good_metadata
    path 'pass_qc.tsv', emit: good_assemblies_list

    script:
    """
    # Create list of assemblies that pass qc & Copy them
    awk -F "\t" '{ if((\$2 < $params.ctg_count) && (\$8 > $params.as_ln_lwr && \$8 < $params.as_ln_upr) && (\$15 > $params.largest_ctg) && (\$17 > $params.gc_lwr && \$17 < $params.gc_upr) && (\$18 > $params.n50)) { print } }' $qcs > pass_qc.tsv
    for assembly_name in `cut -f1 pass_qc.tsv`
    do 
        ln -s ${params.assemblypath}/\${assembly_name}.${fileextension} \${assembly_name}.${fileextension}
    done

    # Create list from good assemblies & create a new metadata file of the filtered results
    ls *.${fileextension} > good_assemblies.txt
    awk -F',' 'NR==FNR{c[\$1]++;next};c[\$1] > 0' good_assemblies.txt $metadata_file > good_metadata.csv
    # Add headers to new metadata file by copying the first line of the orginal metadata file
    printf '%s\n' '0r !head -n 1 $metadata_file' x | ex good_metadata.csv 

    # Remove extensions from metadata file
    awk '{ gsub(/\\.${fileextension}/,"", \$1); print }' good_metadata.csv  > good_noext_metadata.csv
    """
}

// Annotate via prokka. Reference protein file required
process prokka_annotation {
    publishDir  "${params.outdir}/annotation", mode: 'copy', overwrite: true

    input:
    path good_qc_assembly 
    path prokka_ref_file

    output:
    path "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.gff", emit:annotation
    path "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.faa", emit:transl_protein
    path "${good_qc_assembly.baseName}/${good_qc_assembly.baseName}.fna", emit:nucleotide

    script:
    """
    prokka --prefix "${good_qc_assembly.baseName}" --force --proteins ${prokka_ref_file} --centre X --compliant ${good_qc_assembly}
    """
}

// Run Amrfinder using all possible inputs (assembly, ggff, and trans portein files)
process amrfinder {

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

// Run roary as the output directory is needed by piggy as an input
process roary {

    input:
    path annotations //from annotation_roary.collect()

    output:
    path "roary_out/" 

    script:
    """
    roary -v -s -f "roary_out" ${annotations}
    """
}

// Run piggy (extracts intergenic regions)
process piggy {
    publishDir  "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path annotations // from annotation_piggy.collect()
    path roary_dir 

    output:
    path "./piggy_out/IGR_presence_absence.csv", emit: piggy_csv
    path "./piggy_out/IGR_presence_absence.Rtab", emit: piggy_rtab

    script:
    """
    piggy -r "${roary_dir}" -i . -t $params.threads
    """
}

// scoary filtering - make traitfile(Assembly,[Hosts]]) via r script
process gen_scoary_traitfile {
    publishDir  "${params.outdir}/", mode: 'copy', overwrite: true

    input:
    path metadata
    path scoary_datagen_file

    output:
    path "./scoary_traitfile.csv" 

    script:
    """
    Rscript --vanilla $scoary_datagen_file $metadata $params.assembly_column $params.host_column
    """
}

// Filter piggy output (igrs) with scoary
process scoary_igr {
    publishDir  "${params.outdir}/scoary_igr", mode: 'copy', overwrite: true

    input:
    path igr_coll_traitfile
    path igrs 

    output:
    path '*.results.csv' 

    script:
    """
    scoary -t $igr_coll_traitfile -g $igrs --no-time --collapse -p 1.0
    """
}

// Run panaroo 
process panaroo {
    publishDir  "${params.outdir}/panaroo_out", mode: 'copy', overwrite: true

    input:
    path annotations //from annotation_panaroo.collect()

    output:
    path "panaroo_out/gene_presence_absence_roary.csv", emit: pv_csv
    path "panaroo_out/gene_presence_absence.Rtab", emit: pv_rtab

    script:
    """
    panaroo -i ${annotations} -o panaroo_out/ --clean-mode moderate -t $params.threads
    """
}

// Filter panaroo output (pvs) with scoary
process scoary_pv {
    publishDir  "${params.outdir}/scoary_pv", mode: 'copy', overwrite: true

    input:
    path pv_coll_traitfile
    path pvs 

    output:
    path '*.results.csv' 

    script:
    """
    scoary -t $pv_coll_traitfile -g $pvs --no-time --collapse -p 1.0
    """
}

// snippy
process snippy {
    input:
    path assembly 
    path snp_ref_file

    output:
    path "./${assembly.baseName}/" 

    script:
    """
    snippy --prefix "${assembly.baseName}" --outdir "${assembly.baseName}" --ref  "${snp_ref_file}" --ctgs "${assembly}" --force
    """
}

// Get core SNPs
process snippy_core {
    publishDir  "${params.outdir}/snippy_core_out", mode: 'copy', overwrite: true

    input:
    path snippy_folders
    path snp_ref_file

    output:
    path 'core.tab', emit: core_snps
    path 'core.aln', emit: aligned_snps

    script:
    """
    snippy-core --prefix core --ref "${snp_ref_file}" ${snippy_folders}
    """
}

// Get SNP distance of assemblies to use in nonclonal filtering
process snp_dists {

    input:
    path snp_core_aln 

    output:
    path 'snpdist_base.tsv' 

    script:
    """
    snp-dists -m $snp_core_aln > snpdist_base.tsv
    """
}

// Filter for snps
// HACK If no clusters are detected, no list file is produced. so a blank one has been used for now
process clonal_detection {
    publishDir  "${params.outdir}/clonal_detection", mode: 'copy', overwrite: true, pattern : "*.png"

    input:
    path clonal_detection_script
    path metadata
    path snp_dist_file 

    output:
    path "*.list"

    script:
    """
    touch empty.list
    Rscript --vanilla $clonal_detection_script $snp_dist_file $metadata $params.assembly_column $params.host_column $params.year_collected $params.region_column $params.snp_dist_threshold
    """
}

process clonal_filtering {
    publishDir  "${params.outdir}/clonal_detection", mode: 'copy', overwrite: true, pattern : "*.txt"
    
    input:
    path clusters
    path good_assemblies_list
    path metadata

    output:
    path "good_nonclonal_metadata.csv"

    script:
    """
    cut -f1 ${good_assemblies_list} > good_assemblies_all.lst
    cat *.list > cluster_all.txt

    # If no clonal clusters were detected, exit
    if ! [ -s "cluster_all.txt" ];then
        cat $metadata > good_nonclonal_metadata.csv 
        exit
    fi

    # Pick 1 random line from the input file (use as chsen seqs)
    for file in *.list; 
    do 
        head -\$((\${RANDOM} % `wc -l < "\$file"` + 1)) "\$file" | tail -1
    done > cluster_chosen.txt

    # take intersection of above filen
    grep -Fxvf cluster_chosen.txt cluster_all.txt > cluster_discarded.txt
    
    # generate list of all sequences without discarded
    grep -vf cluster_discarded.txt  good_assemblies_all.lst  > all_nonclonal.txt

    # Create a new metadata file of the filtered results & add headers
    awk -F',' 'NR==FNR{c[\$1]++;next};c[\$1] > 0' all_nonclonal.txt $metadata > good_nonclonal_metadata.csv 
    printf '%s\n' '0r !head -n 1 $metadata' x | ex good_nonclonal_metadata.csv 
    """
}

//Concantenate all amr files
process amr_collect {
    input:
    path amr_files

    output:
    path "amr_all.tsv"

    script:
    """
    for file in *.tsv
    do
        name=\${file%_amr.tsv}
        sed -i '1d' \$file
        sed "s/\$/\t\$name/" "\$file"
    done > amr_all.tsv

    # Add headers
    sed -i.bak 1i"Name\tProtein identifier\tContig id\tStart\tStop\tStrand\tGene symbol\tSequence name\tScope\tElement.type\tElement subtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference sequence\t% Identity to reference sequence\tAlignment length\tAccession of closest sequence\tName of closest sequence\tHMM id\tHMM description\tFilename" amr_all.tsv
    """
}

// Process AMR Data for model generation
process amr_process {
    publishDir  "${params.outdir}/model_input", mode: 'copy', overwrite: true, pattern: "*.tsv"

    input:
    path amr_process_script
    path amr_file_all
    path good_nonclonal_metadata
    
    output:
    path "amr_gene_all.tsv", emit: amr_gene_all
    path "amr_gene_bps.tsv", emit: amr_gene_bps
    path "amr_gene_human.tsv", emit: amr_gene_human
    path "amr_class_all.tsv", emit: amr_class_all
    path "amr_class_bps.tsv", emit: amr_class_bps
    path "amr_class_human.tsv", emit: amr_class_human

    script:
    """
    Rscript --vanilla $amr_process_script ${amr_file_all} ${good_nonclonal_metadata} $params.assembly_column $params.host_column
    """
}

//IGR Processing
process igr_process {
    publishDir  "${params.outdir}/model_input", mode: 'copy', overwrite: true, pattern: "*.tsv"

    input:
    path igr_process_script
    path igr_file_all
    path scoary_files
    path good_nonclonal_metadata
    
    output:
    path "igr_all.tsv", emit: igr_all
    path "igr_bps.tsv", emit: igr_bps
    path "igr_human.tsv", emit: igr_human

    script:
    """
    # Combine all scoary files into one & add headers
    for file in *.results.csv
    do
        sed '1d' \$file >>  scoary_all.csv
    done 
    sed -i "1s/^/\$(head -n1 \$file)\\n/" scoary_all.csv

    Rscript --vanilla $igr_process_script ${igr_file_all} scoary_all.csv ${good_nonclonal_metadata} $params.assembly_column $params.host_column
    """
}

// PV Processing
process pv_process {
    publishDir  "${params.outdir}/model_input", mode: 'copy', overwrite: true, pattern: "*.tsv"

    input:
    path pv_process_script
    path pv_file_all
    path scoary_files
    path good_nonclonal_metadata
    
    output:
    path "pv_all.tsv", emit: pv_all
    path "pv_bps.tsv", emit: pv_bps
    path "pv_human.tsv", emit: pv_human

    script:
    """
    # Combine all scoary files into one & add headers
    for file in *.results.csv
    do
        sed '1d' \$file >>  scoary_all.csv
    done 
    sed -i "1s/^/\$(head -n1 \$file)\\n/" scoary_all.csv

    Rscript --vanilla $pv_process_script ${pv_file_all} scoary_all.csv ${good_nonclonal_metadata} $params.assembly_column $params.host_column
    """
}

// SNP processing
process snp_process {
    publishDir  "${params.outdir}/model_input", mode: 'copy', overwrite: true, pattern: "*.tsv"

    input:
    path snps_process_script
    path snp_core
    path good_nonclonal_metadata
    
    output:
    path "snp_abudance_all.tsv", emit: snp_abudance_all
    path "snp_abudance_bps.tsv", emit: snp_abudance_bps
    path "snp_abudance_human.tsv", emit: snp_abudance_human

    script:
    """
    Rscript --vanilla $snps_process_script $snp_core $good_nonclonal_metadata $params.assembly_column $params.host_column
    """
}

// Model generation - Host/Source Attribution
process model_building {
    publishDir "${params.outdir}/models_out/models", mode: 'copy', overwrite: true, pattern: "*.rds"
    publishDir "${params.outdir}/models_out/predictions", mode: 'copy', overwrite: true, pattern: "*.csv"
    publishDir "${params.outdir}/models_out/plots", mode: 'copy', overwrite: true, pattern: "*.png"

    input:
    path amr_class_all
    path amr_class_bps
    path amr_gene_all
    path amr_gene_bps
    path pv_all
    path pv_bps
    path igr_all
    path igr_bps
    path snp_abudance_all
    path snp_abudance_bps
    path model_building_script

    output:
    path "*.rds", emit: models
    path "*.csv", emit: predictions
    path "*.png", emit: plots

    script:
    """
    Rscript --vanilla $model_building_script $amr_class_all $amr_class_bps $amr_gene_all $amr_gene_bps $pv_all $pv_bps $igr_all $igr_bps $snp_abudance_all $snp_abudance_bps $params.model_threads
    """
}

// Model generation - Human scoring models
process model_building_human {
    publishDir "${params.outdir}/models_out/models", mode: 'copy', overwrite: true, pattern: "*.rds"
    publishDir "${params.outdir}/models_out/predictions", mode: 'copy', overwrite: true, pattern: "*.csv"
    publishDir "${params.outdir}/models_out/plots", mode: 'copy', overwrite: true, pattern: "*.png"

    input:
    path amr_class_human
    path amr_gene_human
    path pv_human
    path igr_human
    path snp_abudance_human
    path model_building_human_script

    output:
    path "*.rds", emit: models
    path "*.csv", emit: predictions
    path "*.png", emit: plots

    script:
    """
    Rscript --vanilla $model_building_human_script $amr_class_human $amr_gene_human $pv_human $igr_human $snp_abudance_human $params.model_threads
    """
}

// TODO Add plylogeny (separate workflow)

///////////////////////////////     WORKFLOW    ////////////////////////////////////////////////////////////////////////////////////////////
assemblies = Channel.fromPath("${params.assemblypath}/*.${fileextension}")
metadata_file = file(params.hostdata)

workflow {
    assembly_qc(assemblies)
    printqc(assembly_qc.out.collectFile(keepHeader:true, skip:1), metadata_file)

    prokka_annotation(printqc.out.good_assemblies.flatten(), prokka_ref_file)
    amrfinder(prokka_annotation.out.annotation, prokka_annotation.out.transl_protein, prokka_annotation.out.nucleotide)
    roary(prokka_annotation.out.annotation.collect())
    piggy(prokka_annotation.out.annotation.collect(), roary.out)
    gen_scoary_traitfile(printqc.out.good_metadata, scoary_datagen_file)
    scoary_igr(gen_scoary_traitfile.out, piggy.out.piggy_csv)
    panaroo(prokka_annotation.out.annotation.collect())
    scoary_pv(gen_scoary_traitfile.out, panaroo.out.pv_csv)
    snippy(printqc.out.good_assemblies.flatten(), snp_ref_file)
    snippy_core(snippy.out.collect(), snp_ref_file)

    snp_dists(snippy_core.out.aligned_snps)
    clonal_detection(clonal_detection_script, printqc.out.good_metadata, snp_dists.out)
    clonal_filtering(clonal_detection.out, printqc.out.good_assemblies_list, printqc.out.good_metadata)

    amr_collect(amrfinder.out.collect())
    amr_process(amr_process_script, amr_collect.out, clonal_filtering.out)
    igr_process(igr_process_script, piggy.out.piggy_rtab, scoary_igr.out, clonal_filtering.out)
    pv_process(pv_process_script, panaroo.out.pv_rtab, scoary_pv.out, clonal_filtering.out)
    snp_process(snps_process_script, snippy_core.out.core_snps, clonal_filtering.out)

    model_building( amr_process.out.amr_class_all, amr_process.out.amr_class_bps,
                    amr_process.out.amr_gene_all, amr_process.out.amr_gene_bps,
                    pv_process.out.pv_all, pv_process.out.pv_bps,
                    igr_process.out.igr_all, igr_process.out.igr_bps,
                    snp_process.out.snp_abudance_all, snp_process.out.snp_abudance_bps,
                    model_building_script
                    )
    model_building_human( amr_process.out.amr_class_human,
                    amr_process.out.amr_gene_human, 
                    pv_process.out.pv_human,
                    igr_process.out.igr_human, 
                    snp_process.out.snp_abudance_human, 
                    model_building_human_script
                    )
}

// TODO Break up into subworkflows (one where it doesnt check for clonal)
/*
workflow model_prep {
    if( params.check_clonal )  // defaults to true 
        bar(params.data)
    else
        bar(foo()) }
*/
// TODO Add workflow for inputting new data
