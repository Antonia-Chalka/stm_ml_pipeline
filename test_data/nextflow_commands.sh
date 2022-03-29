# To build models from test data (model_build_in)
nextflow run pipeline_dsl2.nf --assemblypath /home/annita/3.pipeline_data/eb1_run/assemblies/eb1 --hostdata /home/annita/3.pipeline_data/eb1_run/all_metadata.csv --outdir /home/annita/3.pipeline_data/eb1_run/out --year_collected="Collection.Year" --scoarycutoff=0.05 -profile docker -resume

# placeholder for testing defaults
#TODO ADD AND UPDATE CONFIG FILE,  HELP MESSAGE & README
nextflow run pipeline_test.nf 

# Run pipeline on st313 models
nextflow run pipeline_test.nf --assemblypath /home/annita/3.pipeline_data/eb1_run/assemblies/st313 --outdir /home/annita/3.pipeline_data/eb1_run/out/testing_st313 --pv_fasta /home/annita/3.pipeline_data/eb1_run/out/4.model/model_ref/pv_filter.fasta --igr_fasta /home/annita/3.pipeline_data/eb1_run/out/4.model/model_ref/igr_filter.fasta --models /home/annita/3.pipeline_data/eb1_run/out/4.model/models  --snp_core_ref /home/annita/3.pipeline_data/eb1_run/out/2.genomic_features/snp_core_out/core.ref.tab --year_collected="Collection.Year" -profile docker -resume
