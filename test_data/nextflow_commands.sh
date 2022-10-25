# To build models from test data (model_build_in)

# tro test with default dataset


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

