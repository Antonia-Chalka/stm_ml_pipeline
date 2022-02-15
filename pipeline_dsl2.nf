#!/usr/bin/env nextflow
nextflow.enable.dsl=2

////////////////////////////////////////////// HELP SECTION //////////////////////////////////////////////
def helpMessage() {
    log.info """
        Usage:
        The options for running the pipeline to build models are structured as follows:
        --option                       Description/Notes [default value]

    Mandatory Parameters:
        --assemblypath                Directory of your fasta files (full path required) [./test_data/model_build_inut_test]
        --hostdata                    CSV file containing assembly filename, host, year and region [./test_data/model_build_in/metadata.csv]

    Optional Parameters:
        --outdir                       Output directory for models & other data [./out]
    
    Hostdata File Parameters:
        --assembly_column              Column name of your assemblies. Must contain extension eg myassembly.fasta ["Filename"]
        --host_column                  Column name of your hosts ["Source.Host"]
        --region_column                Column of region of origin. Used to detect clonal clusters. Can be empty, but must exist. ["Region"]
        --year_collection              Column of year each assembly was collected. Used to detect clonal clusters. Can be empty, but must exist. ["Year"]

    Assembly Quality Parameters:
        --as_ln_upr                    Maximum accepted assembly length [6000000]
        --as_ln_lwr                    Minumum accepted assembly length [4000000]
        --ctg_count                    Minimum accepted number of contigs [500]
        --largest_ctg                  Minimum accepted length of largest contig [100000] 
        --n50                          Minimum accepted n50 [50000]
        --gc_upr                       Minimum accepted GC % [54]
        --gc_lwr                       Maximum accepted GC % [50]

    Clonal Detection Parameters:
        --snp_dist_threshold           SNP difference used to detect clonal clusters. Used in conjuction with region & collection year metadata. [10]
        
    Tool-specific Parameters:
        --prokka_ref                   'Trusted' protein file for prokka (prokka --proteins) [./data/stm_proteinref.fasta]
        --amr_species                  Assembly species for amrfinder (amrfinder -O) ["Salmonella"]
        --panaroo_mode                 Panaroo assembly filtering mode (panaroo --clean-mode) ["moderate"]
        --snp_ref                      Reference file for snippy (snippy --ref) [./data/stm_sl1344.fasta]
        --threads                      Num of threads to use for panaroo & piggy (panaroo -t & piggy -t) [10]
    
    If you wish to alter the scripts used to generate the models, simply edit the appropriate 'model_building' R scripts in ./data/ - ONLY DO SO IF YOU KNOW WHAT YOU ARE DOING
    """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

////////////////////////////////////////////// MODULES, REFERENCES & SCRIPT FILES //////////////////////////////////////////////
// Modules 
include { assembly_qc                 } from "$projectDir/modules/assembly_qc.nf"
include { printqc                     } from "$projectDir/modules/printqc.nf"
include { prokka_annotation           } from "$projectDir/modules/prokka_annotation.nf"
include { amrfinder                   } from "$projectDir/modules/amrfinder.nf"
include { gen_scoary_traitfile        } from "$projectDir/modules/gen_scoary_traitfile.nf"
include { panaroo                     } from "$projectDir/modules/panaroo.nf"
include { piggy                       } from "$projectDir/modules/piggy.nf"
include { scoary_pv                   } from "$projectDir/modules/scoary_pv.nf"
include { scoary_igr                  } from "$projectDir/modules/scoary_igr.nf"
include { snippy                      } from "$projectDir/modules/snippy.nf"
include { snippy_core                 } from "$projectDir/modules/snippy_core.nf"
include { snp_dists                   } from "$projectDir/modules/snp_dists.nf"
include { clonal_detection            } from "$projectDir/modules/clonal_detection.nf"
include { clonal_filtering            } from "$projectDir/modules/clonal_filtering.nf"
include { amr_collect                 } from "$projectDir/modules/amr_collect.nf"
include { amr_process                 } from "$projectDir/modules/amr_process.nf"
include { igr_process                 } from "$projectDir/modules/igr_process.nf"
include { pv_process                  } from "$projectDir/modules/pv_process.nf"
include { snp_process                 } from "$projectDir/modules/snp_process.nf"
include { model_building              } from "$projectDir/modules/model_building.nf"
include { model_building_human        } from "$projectDir/modules/model_building_human.nf"
include { get_igr_fastas              } from "$projectDir/modules/get_igr_fastas.nf"
include { get_pv_fastas               } from "$projectDir/modules/get_pv_fastas.nf"

// Reference files
prokka_ref_file = file(params.prokka_ref)
snp_ref_file = file(params.snp_ref)

// R script files
scoary_datagen_file=file("$projectDir/data/scoary_generate_tabfile.R")
clonal_detection_script=file("$projectDir/data/filter_script.R")
amr_process_script=file("$projectDir/data/input_amr.R")
pv_process_script=file("$projectDir/data/input_pv.R")
igr_process_script=file("$projectDir/data/input_igr.R")
snps_process_script=file("$projectDir/data/input_snps.R")
model_building_script=file("$projectDir/data/model_building.R")
model_building_human_script=file("$projectDir/data/model_building_human.R")

//////////////////////////////////////////////     WORKFLOW    //////////////////////////////////////////////
assemblies = Channel.fromPath("${params.assemblypath}/*.${params.fileextension}")
metadata_file = file(params.hostdata)

workflow {
    assembly_qc(assemblies)
    printqc(
        assembly_qc.out.collectFile(keepHeader:true, skip:1, sort:true,storeDir:"${params.outdir}/1.assembly_quality"), 
        metadata_file)
    gen_scoary_traitfile(
        printqc.out.good_metadata, 
        scoary_datagen_file)

    prokka_annotation(
        printqc.out.good_assemblies.flatten(), 
        prokka_ref_file)
    amrfinder(
        prokka_annotation.out.annotation, 
        prokka_annotation.out.transl_protein,
        prokka_annotation.out.nucleotide)

    panaroo(prokka_annotation.out.annotation.collect())
    piggy(
        prokka_annotation.out.annotation.collect(), 
        panaroo.out.roary_dir)
    scoary_pv(
        gen_scoary_traitfile.out, 
        anaroo.out.pv_csv)
    scoary_igr(
        gen_scoary_traitfile.out, 
        piggy.out.piggy_csv)

    snippy(
        printqc.out.good_assemblies.flatten(), 
        snp_ref_file)
    snippy_core(
        snippy.out.collect(), 
        snp_ref_file)
    snp_dists(snippy_core.out.aligned_snps)
    clonal_detection(
        clonal_detection_script, 
        printqc.out.good_metadata, 
        snp_dists.out)
    clonal_filtering(
        clonal_detection.out.clusters, 
        printqc.out.good_assemblies_list, 
        printqc.out.good_metadata)

    amr_collect(amrfinder.out.collect())
    amr_process(
        amr_process_script, 
        amr_collect.out, 
        clonal_filtering.out)
    igr_process(
        igr_process_script, 
        piggy.out.piggy_rtab, 
        scoary_igr.out, 
        clonal_filtering.out)
    pv_process(
        pv_process_script, 
        panaroo.out.pv_rtab, 
        scoary_pv.out, 
        clonal_filtering.out)
    snp_process(
        snps_process_script, 
        snippy_core.out.core_snps,
        clonal_filtering.out)

    get_igr_fastas(
        igr_process.out.igr_all, 
        piggy.out.piggy_ref)
    get_pv_fastas(
        pv_process.out.pv_all,
        panaroo.out.pv_ref)

    model_building( 
        amr_process.out.amr_class_all, 
        amr_process.out.amr_class_bps,
        amr_process.out.amr_gene_all, 
        amr_process.out.amr_gene_bps,
        pv_process.out.pv_all, 
        pv_process.out.pv_bps,
        igr_process.out.igr_all, 
        igr_process.out.igr_bps,
        snp_process.out.snp_abudance_all, 
        snp_process.out.snp_abudance_bps,
        model_building_script
                    )
    model_building_human( 
        amr_process.out.amr_class_human,
        amr_process.out.amr_gene_human, 
        pv_process.out.pv_human,
        igr_process.out.igr_human, 
        snp_process.out.snp_abudance_human, 
        model_building_human_script
        )

    
}
