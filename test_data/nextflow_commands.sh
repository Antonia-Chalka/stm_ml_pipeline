# To build models from test data (model_build_in)
nextflow run pipeline_dsl2_subworkflow.nf \
    --assemblypath /home/annita/3.pipeline_data/stm_ml_pipeline-1/test_data/model_build_in \
    --hostdata /home/annita/3.pipeline_data/stm_ml_pipeline-1/test_data/model_build_in/metadata.csv \
    --outdir /home/annita/3.pipeline_data/stm_ml_pipeline-1/test_data/model_build_out \
    --year_collected="Year" \
    --scoarycutoff=0.05 \
    --prokka_extra "--genus Enterococcus" \
    -profile docker \
    -resume

# Use test data dn only get genomic features that do not require you have the entire dataset
nextflow run pipeline_dsl2_subworkflow.nf \
    --assemblypath /home/annita/3.pipeline_data/stm_ml_pipeline-1/test_data/model_build_in \
    --hostdata /home/annita/3.pipeline_data/stm_ml_pipeline-1/test_data/model_build_in/metadata.csv \
    --outdir /home/annita/3.pipeline_data/stm_ml_pipeline-1/test_data/model_build_out \
    --year_collected="Year" \
    --scoarycutoff=0.05 \
    --prokka_extra "--genus Enterococcus" \
    -profile docker \
    -entry incomplete_dataset \
    -resume

    # June analysis
nextflow run pipeline_dsl2_subworkflow.nf \
    --assemblypath /mnt/i/june_analysis/assemblies \
    --hostdata /mnt/i/june_analysis/additional_june.csv \
    --outdir //mnt/i/june_analysis/assemblies/out \
    --assembly_column "Filename" \
    --host_column "Source Type" \
    --year_collection "Collection Year" \
    --scoarycutoff=0.05 \
    --prokka_extra "--genus Enterococcus" \
    -profile docker \
    -entry incomplete_dataset \
    -resume

#
nextflow run pipeline_dsl2_subworkflow.nf \
    --assemblypath /mnt/i/final_analysis/1.assemblies_in \
    --prokka_ref /mnt/i/final_analysis/ecoli_trusted_uniprot.fasta  \   
    --snp_ref /mnt/i/final_analysis/EC958_SNIPPYREF.chr.fa   \  
    --hostdata /mnt/i/final_analysis/1.ecoli_met.csv     \
    --outdir /mnt/i/final_analysis/2.pipeline_out     \
    --assembly_column "Filename"     \
    --host_column "Host"     \
    --year_collected "Collection Year"     \
    --scoarycutoff=0.05     \
    --prokka_extra "--genus Enterococcus"     \
    -profile docker     \
    -entry feature_create     \
    -resume



# tro test with default dataset
nextflow run pipeline_dsl2_subworkflow.nf \
    --assemblypath /mnt/i/final_analysis/1.assemblies_in \
    --prokka_ref /mnt/i/final_analysis/ecoli_trusted_uniprot.fasta \
    --snp_ref /mnt/i/final_analysis/EC958_SNIPPYREF.chr.fa \
    --hostdata /mnt/i/final_analysis/1.ecoli_met.csv \
    --outdir /mnt/i/final_analysis/2.pipeline_out \
    --assembly_column "Filename" \
    --host_column "Host" \
    --year_collected "Collection Year" \
    --scoarycutoff=0.05 \
    --prokka_extra "--genus Enterococcus" \
    -profile docker \
    -entry feature_create \
    -resume
    



# EB1 Build
nextflow run pipeline_dsl2.nf --assemblypath /home/annita/3.pipeline_data/eb1_run/assemblies/eb1 --hostdata /home/annita/3.pipeline_data/eb1_run/all_metadata.csv --outdir /home/annita/3.pipeline_data/eb1_run/out --year_collected="Collection.Year" --scoarycutoff=0.05 --prokka_extra "--genus Enterococcus" -profile docker -resume

# Run pipeline on other datasets
nextflow run pipeline_test.nf --assemblypath /home/annita/3.pipeline_data/eb1_run/assemblies/st313 --outdir /home/annita/3.pipeline_data/eb1_run/out/testing_st313 --pv_fasta /home/annita/3.pipeline_data/eb1_run/out/4.model/model_ref/pv_filter.fasta --igr_fasta /home/annita/3.pipeline_data/eb1_run/out/4.model/model_ref/igr_filter.fasta --models /home/annita/3.pipeline_data/eb1_run/out/4.model/models  --snp_core_ref /home/annita/3.pipeline_data/eb1_run/out/2.genomic_features/snp_core_out/core.ref.tab --year_collected="Collection.Year" -profile docker -resume

nextflow run pipeline_test.nf --assemblypath /home/annita/3.pipeline_data/eb1_run/assemblies/outbreak_t5.7665 --outdir /home/annita/3.pipeline_data/eb1_run/out/outbreak_t5.7665 --pv_fasta /home/annita/3.pipeline_data/eb1_run/out/4.model/model_ref/pv_filter.fasta --igr_fasta /home/annita/3.pipeline_data/eb1_run/out/4.model/model_ref/igr_filter.fasta --models /home/annita/3.pipeline_data/eb1_run/out/4.model/models  --snp_core_ref /home/annita/3.pipeline_data/eb1_run/out/2.genomic_features/snp_core_out/core.ref.tab --year_collected="Collection.Year" -profile docker -resume


# New USA data (2020-2022)
nextflow run pipeline_test.nf --assemblypath /mnt/f/5.additional_datasets/stm_usa_2020-2022/usa-20-22-bps/in \
    --outdir /mnt/f/3.pipeline_data/eb1_run/out/usa_new_20-22 \
    --pv_fasta /mnt/f/3.pipeline_data/eb1_run/out_old/4.model/model_ref/pv_filter.fasta \
    --igr_fasta /mnt/f/3.pipeline_data/eb1_run/out_old/4.model/model_ref/igr_filter.fasta \
    --models /mnt/f/3.pipeline_data/eb1_run/out_old/4.model/models \
    --snp_core_ref /mnt/f/3.pipeline_data/eb1_run/out_old/2.genomic_features/snp_core_out/core.ref.tab \
    --assembly_column="Filename" \
    -profile docker 


# phage run
nextflow run pipeline_dsl2.nf --assemblypath /mnt/f/Phages/models/phage_interaction_cleaned \
    --hostdata /mnt/f/Phages/models/test.csv \
    --prokka_ref /mnt/f/Phages/models/ecoli_trusted_uniprot.fasta \
    --snp_ref /mnt/f/Phages/models/EC958.chr.fa \
    --outdir /mnt/f/Phages/models/test_out \
    --scoary_threshold 70 \
    --scoary_gen_file ~/3.pipeline_data/stm_ml_pipeline-1/data/scoary_datagen_phage.R \
    --scoarycutoff=0.05 \
    --prokka_extra "--genus Escherichia" \
    --assembly_column  "Isolate" \
    --host_column "CHAP1" \
    --as_ln_upr 6000000 \
    --as_ln_lwr 4000000 \
    --ctg_count 600 \
    --largest_ctg 1 \
    --n50 1 \
    --gc_upr 99 \
    --gc_lwr 1 \
    -profile docker \
    -resume

#phage test



# vesa run
nextflow run pipeline_dsl2.nf --assemblypath /mnt/f/Other_Projects/Vesa_ecoli/model_in/fasta_in \
    --hostdata /mnt/f/Other_Projects/Vesa_ecoli/model_in/urban_zoo_metadata_datprepped_clean.csv \
    --prokka_ref /mnt/f/Other_Projects/Vesa_ecoli/model_in/ecoli_trusted_uniprot.fasta \
    --snp_ref /mnt/f/Other_Projects/Vesa_ecoli/model_in/EC958.chr.fa \
    --outdir /mnt/f/Other_Projects/Vesa_ecoli/model_out \
    --scoarycutoff=0.05 \
    --assembly_column="sample_name" \
    --host_column="Source_Level1" \
    --region_column="Sublocation" \
    --year_collection="Year" \
    --amr_species="Escherichia" \
    --prokka_extra "--genus Escherichia" \
    --as_ln_upr 6000000 \
    --as_ln_lwr 4000000 \
    --ctg_count 500 \
    --largest_ctg 1 \
    --n50 1 \
    --gc_upr 99 \
    --gc_lwr 1 \
    -profile docker \
    -resume



