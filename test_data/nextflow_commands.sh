# To build models from test data (model_build_in)

nextflow run pipeline_dsl2.nf --outdir /home/annita/repos/stm_ml_pipeline-1/test_data/model_build_out/ --assemblypath /home/annita/repos/stm_ml_pipeline-1/test_data/model_build_in --hostdata /home/annita/repos/stm_ml_pipeline-1/test_data/model_build_in/metadata.csv --year_collected="Collection.Year" -profile docker -resume

nextflow run pipeline_dsl2.nf -profile docker -resume

nextflow run pipeline_dsl2.nf --outdir /home/annita/repos/initial_data_run/out/ --assemblypath /home/annita/repos/initial_data_run/assemblies --hostdata /home/annita/repos/initial_data_run/all_metadata.csv --year_collected="Collection.Year" -profile docker -resume


# Run pipeline that tests human seq against old models
nextflow run pipeline_test.nf --assemblypath /home/annita/repos/initial_data_run/old_model_pred/in/human_assemblies --hostdata /home/annita/repos/initial_data_run/old_model_pred/in/all_metadata.csv --outdir /home/annita/repos/initial_data_run/old_model_pred/out --pv_fasta /home/annita/repos/initial_data_run/old_model_pred/in/pv_filter.fasta --igr_fasta /home/annita/repos/initial_data_run/old_model_pred/in/igr_filter.fasta --models /home/annita/repos/initial_data_run/old_model_pred/in/models/ --year_collected="Collection.Year" -profile docker -resume
