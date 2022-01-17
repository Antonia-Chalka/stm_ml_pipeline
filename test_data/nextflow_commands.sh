# To build models from test data (model_build_in)

nextflow run pipeline_dsl2.nf --outdir /home/annita/repos/stm_ml_pipeline-1/test_data/model_build_out/ --assemblypath /home/annita/repos/stm_ml_pipeline-1/test_data/model_build_in --hostdata /home/annita/repos/stm_ml_pipeline-1/test_data/model_build_in/metadata.csv --year_collected="Collection.Year" -profile docker -resume

nextflow run pipeline_dsl2.nf -profile docker -resume
