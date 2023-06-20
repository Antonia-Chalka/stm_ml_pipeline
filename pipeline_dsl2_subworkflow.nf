#!/usr/bin/env nextflow
nextflow.enable.dsl=2

////////////////////////////////////////////// HELP SECTION //////////////////////////////////////////////
def helpMessage() {
    log.info """
        Usage:
        The options for running the pipeline to build models are structured as follows:
        --option                       Description/Notes [default value]

    Mandatory Parameters:
        --assemblypath                  Directory of your fasta files (full path required) [./test_data/model_build_inut_test]
        --hostdata                      CSV file containing assembly filename, host, year and region [./test_data/model_build_in/metadata.csv]

    Optional Parameters:
        --outdir                        Output directory for models & other data [./out]
    
    Hostdata File Parameters:
        --assembly_column               Column name of your assemblies. Must contain extension eg myassembly.fasta ["Filename"]
        --host_column                   Column name of your hosts ["Source.Host"]
        --region_column                 Column of region of origin. Used to detect clonal clusters. Can be empty, but must exist. ["Region"]
        --year_collection               Column of year each assembly was collected. Used to detect clonal clusters. Can be empty, but must exist. ["Year"]

    Assembly Quality Parameters:
        --as_ln_upr                     Maximum accepted assembly length [6000000]
        --as_ln_lwr                     Minimum accepted assembly length [4000000]
        --ctg_count                     Minimum accepted number of contigs [500]
        --largest_ctg                   Minimum accepted length of largest contig [100000] 
        --n50                           Minimum accepted n50 [50000]
        --gc_upr                        Maximum accepted GC % [54]
        --gc_lwr                        Minimum accepted GC % [50]

    Clonal Detection Parameters:
        --snp_dist_threshold            SNP difference used to detect clonal clusters. Used in conjunction with region & collection year metadata. [10]
        
    Tool-specific Parameters:
        --prokka_ref                    'Trusted' protein file for prokka (prokka --proteins) [./data/stm_proteinref.fasta]
        --prokka_extra                  Additional options for prokka anootation, eg."--genus Enterococcus" []
        --amr_species                   Assembly species for amrfinder (amrfinder -O) ["Salmonella"]
        --panaroo_mode                  Panaroo assembly filtering mode (panaroo --clean-mode) ["moderate"]
        --snp_ref                       Reference file for snippy (snippy --ref) [./data/stm_sl1344.fasta]
        --threads                       Num of threads to use for panaroo & piggy (panaroo -t & piggy -t) [10]
    
    Feature Filtering Parameters:
        snp_lower                       Lower threshold abundance of SNPs to be excluded from odel training [0.1]
        snp_upper                       Lower threshold abundance of SNPs to be excluded from odel training [99]
        scoarycutoff                    Threshold for filtering PV/IGRs based on scoary's bonferroni corrected p value [1.1]

    Nextflow Options
        profile                         Run with docker or singularity (options are docker or singularity)

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
include { assembly_qc                                   } from "$projectDir/modules/assembly_qc.nf"
include { printqc                                       } from "$projectDir/modules/printqc.nf"
include { prokka_annotation                             } from "$projectDir/modules/prokka_annotation.nf"
include { amrfinder                                     } from "$projectDir/modules/amrfinder.nf"
include { gen_scoary_traitfile                          } from "$projectDir/modules/gen_scoary_traitfile.nf"
include { panaroo                                       } from "$projectDir/modules/panaroo.nf"
include { piggy                                         } from "$projectDir/modules/piggy.nf"
include { scoary_pv                                     } from "$projectDir/modules/scoary_pv.nf"
include { scoary_igr                                    } from "$projectDir/modules/scoary_igr.nf"
include { snippy                                        } from "$projectDir/modules/snippy.nf"
include { snippy_core                                   } from "$projectDir/modules/snippy_core.nf"
include { snp_dists                                     } from "$projectDir/modules/snp_dists.nf"
include { clonal_detection                              } from "$projectDir/modules/clonal_detection.nf"
include { clonal_filtering                              } from "$projectDir/modules/clonal_filtering.nf"
include { amr_collect                                   } from "$projectDir/modules/amr_collect.nf"
include { amr_process                                   } from "$projectDir/modules/amr_process.nf"
include { igr_process                                   } from "$projectDir/modules/igr_process.nf"
include { pv_process                                    } from "$projectDir/modules/pv_process.nf"
include { snp_process                                   } from "$projectDir/modules/snp_process.nf"
include { model_building as model_building_amr          } from "$projectDir/modules/model_building.nf"
include { model_building as model_building_pv           } from "$projectDir/modules/model_building.nf"
include { model_building as model_building_igr          } from "$projectDir/modules/model_building.nf"
include { model_building as model_building_snp          } from "$projectDir/modules/model_building.nf"
include { get_igr_fastas                                } from "$projectDir/modules/get_igr_fastas.nf"
include { get_pv_fastas                                 } from "$projectDir/modules/get_pv_fastas.nf"

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
model_building_script=file("$projectDir/data/single_model_build.R")

//////////////////////////////////////////////     WORKFLOW    //////////////////////////////////////////////
assemblies = Channel.fromPath("${params.assemblypath}/*.${params.fileextension}")
metadata_file = file(params.hostdata)

//////////////////////////


/// update testing ppipeline (??)
// predict which host from human
//[redict if a livestock one is pathogenic

/// SNP with Scoary
/// SNP address
/// Input phylotree to scoary
/// generate phylotree
/// iterative feature reduction input
// r TO SCIKIT
/// OUTPUT VISUALISATIONS

workflow single_seq_operations {
    // filtering, annotation, snp, AMR
    take: 
        assemblies
        metadata_file
    main:
        assembly_qc(assemblies) //ok no ref
        printqc(
            assembly_qc.out.collectFile(keepHeader:true, skip:1, sort:true,storeDir:"${params.outdir}/1.assembly_quality"), 
            metadata_file)
        
        prokka_annotation(
            printqc.out.good_assemblies.flatten(), 
            prokka_ref_file)
        amrfinder(
            prokka_annotation.out.annotation, 
            prokka_annotation.out.transl_protein,
            prokka_annotation.out.nucleotide)
        snippy(
            printqc.out.good_assemblies.flatten(), 
            snp_ref_file)
    emit:
        good_assemblies_list = printqc.out.good_assemblies_list 
        good_metadata = printqc.out.good_metadata
        annnotation = prokka_annotation.out.annotation.collect()
        amrfinder = amrfinder.out.collect()
        snippy = snippy.out.collect()
}

workflow multi_seq_operations {
    // pangenome intergenic core snps, snps filtering
    take: 
        good_assemblies_list // printqc.out.good_assemblies_list 
        good_metadata //printqc.out.good_metadata
<<<<<<< HEAD
        annotation// TODO how to make it take a dir as an aletrnative? //prokka_annotation.out.annotation.collect()
        snippy //snippy.out.collect()
    main:
        gen_scoary_traitfile(
            good_metadata, // change to input metadata file
            scoary_datagen_file)

        panaroo(annotation)  // from prokka
        piggy(
            annotation,  //from prokka
            panaroo.out.roary_dir)

        snippy_core(
        snippy, 
=======
        annnotation// TODO how to make it take a dir as an aletrnative? //prokka_annotation.out.annotation.collect()
        snippy //snippy.out.collect()
    main:
        gen_scoary_traitfile(
            printqc.out.good_metadata, // change to input metadata file
            scoary_datagen_file)

        panaroo(prokka_annotation.out.annotation.collect())  // from prokka
        piggy(
            prokka_annotation.out.annotation.collect(),  //from prokka
            panaroo.out.roary_dir)

        snippy_core(
        snippy.out.collect(), 
>>>>>>> 34a47416a03aa71a06c48d34dc048457c0fee6b9
        snp_ref_file)
        snp_dists(snippy_core.out.aligned_snps)
        clonal_detection(
            clonal_detection_script, 
<<<<<<< HEAD
            good_metadata, 
            snp_dists.out)
        clonal_filtering(
            clonal_detection.out.clusters, 
            good_assemblies_list, // change to input custom good assemblies file (?)
            good_metadata) //change to input metadat file
    emit:
        scoary_traitfile = gen_scoary_traitfile.out
=======
            printqc.out.good_metadata, 
            snp_dists.out)
        clonal_filtering(
            clonal_detection.out.clusters, 
            printqc.out.good_assemblies_list, // change to input custom good assemblies file (?)
            printqc.out.good_metadata) //change to input metadat file
    emit:
>>>>>>> 34a47416a03aa71a06c48d34dc048457c0fee6b9
        pv_csv = panaroo.out.pv_csv
        pv_rtab = panaroo.out.pv_rtab
        pv_ref = panaroo.out.pv_ref
        piggy_csv = piggy.out.piggy_csv
        piggy_rtab = piggy.out.piggy_rtab
        piggy_ref = piggy.out.piggy_ref
        core_snps = snippy_core.out.core_snps
        clonal_filtering_out = clonal_filtering.out
}

workflow phylogeny {
    // snps, accessory genome, 
    // hierbaps gubbins 
    take: data
    main:
        foo(data)
        baz(foo.out)
    emit:
        baz.out
}

workflow filtering {
    take:
<<<<<<< HEAD
        scoary_traitfile // gen_scoary_traitfile.out
=======
>>>>>>> 34a47416a03aa71a06c48d34dc048457c0fee6b9
        amrfinder // amrfinder.out.collect()
        pv_csv //panaroo.out.pv_csv
        pv_rtab //panaroo.out.pv_rtab
        pv_ref //panaroo.out.pv_ref
        piggy_csv //piggy.out.piggy_csv
<<<<<<< HEAD
        piggy_rtab // piggy.out.piggy_rtab
=======
>>>>>>> 34a47416a03aa71a06c48d34dc048457c0fee6b9
        piggy_ref //piggy.out.piggy_ref
        core_snps //snippy_core.out.core_snps
        clonal_filtering_out //clonal_filtering.out
    main:
    // with/without tree
    // pv igr snp
        scoary_pv(
<<<<<<< HEAD
            scoary_traitfile, 
            pv_csv)
        scoary_igr(
            scoary_traitfile, 
            piggy_csv)

        amr_collect(amrfinder)
        amr_process(
            amr_process_script, 
            amr_collect.out, 
            clonal_filtering_out)
        igr_process(
            igr_process_script, 
            piggy_rtab, 
            scoary_igr.out, 
            clonal_filtering_out)
        pv_process(
            pv_process_script, 
            pv_rtab, 
            scoary_pv.out, 
            clonal_filtering_out)
        snp_process(
            snps_process_script, 
            core_snps,
            clonal_filtering_out)

        get_igr_fastas(
            igr_process.out.igr_all, 
            piggy_ref)
        get_pv_fastas(
            pv_process.out.pv_all,
            pv_ref)
=======
            gen_scoary_traitfile.out, 
            panaroo.out.pv_csv)
        scoary_igr(
            gen_scoary_traitfile.out, 
            piggy.out.piggy_csv)

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
>>>>>>> 34a47416a03aa71a06c48d34dc048457c0fee6b9

    emit:
        amr_process_out = amr_process.out.amr_inputs.flatten()
        pv_process_out = pv_process.out.pv_inputs.flatten()
        igr_process_out = igr_process.out.igr_inputs.flatten()
        snp_process_out = snp_process.out.snp_inputs.flatten()
}

workflow model_training {
    // baps fold
    take:
        amr_process_out //amr_process.out.amr_inputs.flatten()
        pv_process_out //pv_process.out.pv_inputs.flatten()
        igr_process_out //igr_process.out.igr_inputs.flatten()
        snp_process_out //snp_process.out.snp_inputs.flatten()

    main:
    // testing training validation
    // scikot
    // shap values?
    model_building_amr( 
        model_building_script,
<<<<<<< HEAD
        amr_process_out)
    model_building_pv( 
        model_building_script,
        pv_process_out)
    model_building_igr( 
        model_building_script,
        igr_process_out)
    model_building_snp( 
        model_building_script,
        snp_process_out)
=======
        amr_process.out.amr_inputs.flatten())
    model_building_pv( 
        model_building_script,
        pv_process.out.pv_inputs.flatten())
    model_building_igr( 
        model_building_script,
        igr_process.out.igr_inputs.flatten())
    model_building_snp( 
        model_building_script,
        snp_process.out.snp_inputs.flatten())
>>>>>>> 34a47416a03aa71a06c48d34dc048457c0fee6b9
}

workflow model_testing {
    // 
}

<<<<<<< HEAD

workflow incomplete_dataset {
    single_seq_operations(assemblies, metadata_file)
}

workflow incomplete_dataset_phylogeny { // TODO add an if statement and merge with the above?
    single_seq_operations(assemblies, metadata_file)
}

workflow {

    single_seq_operations(assemblies, metadata_file)
    multi_seq_operations(
        single_seq_operations.out.good_assemblies_list, //may need to change? wht happns if ffragemnted runs are run through?
        single_seq_operations.out.good_metadata,
        single_seq_operations.out.annnotation,
        single_seq_operations.out.snippy
        )
    filtering(
        multi_seq_operations.out.scoary_traitfile,
        single_seq_operations.out.amrfinder,
        multi_seq_operations.out.pv_csv,
        multi_seq_operations.out.pv_rtab,
        multi_seq_operations.out.pv_ref,
        multi_seq_operations.out.piggy_csv,
        multi_seq_operations.out.piggy_rtab,
        multi_seq_operations.out.piggy_ref,
        multi_seq_operations.out.core_snps,
        multi_seq_operations.out.clonal_filtering_out
=======
workflow {
    single_seq_operations(assemblies, metadata_file)
    multi_seq_operations(
        single_seq_operations.good_assemblies_list, //may need to change? wht happns if ffragemnted runs are run through?
        single_seq_operations.good_metadata,
        single_seq_operations.annnotation,
        single_seq_operations.snippy
        )
    filtering(
       multi_seq_operations.amrfinder 
       multi_seq_operations.pv_csv 
       multi_seq_operations.pv_rtab
       multi_seq_operations.pv_ref
       multi_seq_operations.piggy_csv
       multi_seq_operations.piggy_ref
       multi_seq_operations.core_snps
       multi_seq_operations.clonal_filtering_out
>>>>>>> 34a47416a03aa71a06c48d34dc048457c0fee6b9
    )

    // TODO if statement whether training or testing (default= "train" )
    model_training(
<<<<<<< HEAD
        filtering.out.amr_process_out,
        filtering.out.pv_process_out,
        filtering.out.igr_process_out,
        filtering.out.snp_process_out
=======
        filtering.amr_process_out
        filtering.pv_process_out
        filtering.igr_process_out
        filtering.snp_process_out
>>>>>>> 34a47416a03aa71a06c48d34dc048457c0fee6b9
    )
}


// OLD PIPELINE

// assembly_qc(assemblies)
    // printqc(
    //     assembly_qc.out.collectFile(keepHeader:true, skip:1, sort:true,storeDir:"${params.outdir}/1.assembly_quality"), 
    //     metadata_file)
    // gen_scoary_traitfile(
    //     printqc.out.good_metadata, 
    //     scoary_datagen_file)

    // prokka_annotation(
    //     printqc.out.good_assemblies.flatten(), 
    //     prokka_ref_file)
    // amrfinder(
    //     prokka_annotation.out.annotation, 
    //     prokka_annotation.out.transl_protein,
    //     prokka_annotation.out.nucleotide)

    // panaroo(prokka_annotation.out.annotation.collect())
    // piggy(
    //     prokka_annotation.out.annotation.collect(), 
    //     panaroo.out.roary_dir)
    // scoary_pv(
    //     gen_scoary_traitfile.out, 
    //     panaroo.out.pv_csv)
    // scoary_igr(
    //     gen_scoary_traitfile.out, 
    //     piggy.out.piggy_csv)

    //snippy(
    //    printqc.out.good_assemblies.flatten(), 
    //    snp_ref_file)
    // snippy_core(
    //     snippy.out.collect(), 
    //     snp_ref_file)
    // snp_dists(snippy_core.out.aligned_snps)
    // clonal_detection(
    //     clonal_detection_script, 
    //     printqc.out.good_metadata, 
    //     snp_dists.out)
    // clonal_filtering(
    //     clonal_detection.out.clusters, 
    //     printqc.out.good_assemblies_list, 
    //     printqc.out.good_metadata)

    // amr_collect(amrfinder.out.collect())
    // amr_process(
    //     amr_process_script, 
    //     amr_collect.out, 
    //     clonal_filtering.out)
    // igr_process(
    //     igr_process_script, 
    //     piggy.out.piggy_rtab, 
    //     scoary_igr.out, 
    //     clonal_filtering.out)
    // pv_process(
    //     pv_process_script, 
    //     panaroo.out.pv_rtab, 
    //     scoary_pv.out, 
    //     clonal_filtering.out)
    // snp_process(
    //     snps_process_script, 
    //     snippy_core.out.core_snps,
    //     clonal_filtering.out)

    // get_igr_fastas(
    //     igr_process.out.igr_all, 
    //     piggy.out.piggy_ref)
    // get_pv_fastas(
    //     pv_process.out.pv_all,
    //     panaroo.out.pv_ref)

    // model_building_amr( 
    //     model_building_script,
    //     amr_process.out.amr_inputs.flatten())
    // model_building_pv( 
    //     model_building_script,
    //     pv_process.out.pv_inputs.flatten())
    // model_building_igr( 
    //     model_building_script,
    //     igr_process.out.igr_inputs.flatten())
    // model_building_snp( 
    //     model_building_script,
    //     snp_process.out.snp_inputs.flatten())

