# To build models from test data (model_build_in)

nextflow run pipeline_dsl2.nf --outdir /home/annita/repos/stm_ml_pipeline-1/test_data/model_build_out/ --assemblypath /home/annita/repos/stm_ml_pipeline-1/test_data/model_build_in --hostdata /home/annita/repos/stm_ml_pipeline-1/test_data/model_build_in/metadata.csv --year_collected="Collection.Year" -profile docker -resume

nextflow run pipeline_dsl2.nf -profile docker -resume

nextflow run pipeline_dsl2.nf --outdir /home/annita/repos/initial_data_run/out/ --assemblypath /home/annita/repos/initial_data_run/assemblies --hostdata /home/annita/repos/initial_data_run/all_metadata.csv --year_collected="Collection.Year" -profile docker -resume



sed '/>/ s/_+.*//g' representative_clusters_merged.fasta > cleaned.fasta
head -n 1 igr_all.tsv | sed 's/\"//g' | tr -s ' \t' '\n' | head -n -1 > igr.list
seqtk subseq cleaned.fasta igr.list > igr_filter.fasta


combine assemblies and give unique headers

makeblastdb -in all_assemblies.fasta (file-name of the contigs or whatever) -dbtype nucl  -out assemblyblastdb

tblastx -query input (with the gene file) -db assemblyblastdb -out outname (file name with the results)

