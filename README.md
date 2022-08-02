# Workflow for Host Attribution Machine Learning Models 

A DSL2 Nextflow & Docker pipeline used to build bacterial source attribution machine learning models from assembled genomes, using SNPs, protein variants (PVs), intergenic regions (IGRs) and AMR profiles. Currently, it has only been tested for *Salmonella typhimurium* sequences.

As used in:

TODO Add citation

## Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Building Models](#building-models)
  - [Required Inputs (Model Building)](#required-inputs-model-building)
  - [Advanced Runs (Model Building)](#advanced-runs-model-building)
  - [Full Parameters (Model Building)](#full-parameters-model-building)
  - [Outputs (Model Building)](#outputs-model-building)
  - [Workflow (Model Building)](#workflow-model-building)
  - [Clonal Filtering](#clonal-filtering)
- [Testing Models](#testing-models)
  - [Advanced Runs (Testing Models)](#advanced-runs-testing-models)
  - [Full Parameters (Testing Models)](#full-parameters-testing-models)
  - [Outputs (Testing Models)](#outputs-testing-models)
  - [Workflow (Testing Models)](#workflow-testing-models)
- [Benchmarks](#benchmarks)
- [Citation](#citation)

## Installation

- [Install Java](https://java.com/en/download/help/download_options.html)
- [Install Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- [Install Docker](https://docs.docker.com/get-docker/)
- Pull Required Docker/Singularity Images
  - With Docker:
     ``` bash
      docker pull staphb/quast:5.0.2 ;
      docker pull staphb/prokka:1.14.5 ;
      docker pull staphb/ncbi-amrfinderplus:3.10.5 ;
      docker pull staphb/piggy:1.5 ;
      docker pull quay.io/biocontainers/panaroo:1.2.9--pyhdfd78af_0 ;
      docker pull staphb/snippy:4.6.0 ;
      docker pull rocker/tidyverse:4.0.5 ;
      docker pull quay.io/biocontainers/scoary:1.6.16--py_2 ;
      docker pull staphb/snp-dists:0.8.2 ;
      docker pull annitachalka/r_model_build:1.01 ;
      docker pull ncbi/blast:2.12.0 ;
      docker pull staphb/seqtk:1.3 
      ```
  - With Singularity:
    ``` bash
    singularity pull docker://staphb/quast:5.0.2  ;
    singularity pull docker://staphb/prokka:1.14.5 ; 
    singularity pull docker://staphb/ncbi-amrfinderplus:3.10.5 ; 
    singularity pull docker://staphb/piggy:1.5 ; 
    singularity pull docker://quay.io/biocontainers/panaroo:1.2.9--pyhdfd78af_0 ; 
    singularity pull docker://staphb/snippy:4.6.0 ; 
    singularity pull docker://rocker/tidyverse:4.0.5 ; 
    singularity pull docker://quay.io/biocontainers/scoary:1.6.16--py_2 ; 
    singularity pull docker://staphb/snp-dists:0.8.2 ; 
    singularity pull docker://annitachalka/r_model_build:1.01
    singularity pull docker://ncbi/blast:2.12.0
    singularity pull docker://ncbi/staphb/seqtk:1.3
    ```

## Quick Start

Build models from a set of assemblies:

``` bash
nextflow run pipeline_dsl2.nf --assemblypath "/home/username/myproject/input_data" --hostdata "/home/username/myproject/input_data/metadata.csv" --resume
```

Get predictions of existing models from a set of assemblies:

TODO CHANGE TO ADD MODEL ARGUMENT

``` bash
nextflow run pipeline_dsl2.nf --assemblypath "/home/username/myproject/input_data" --hostdata "/home/username/myproject/input_data/metadata.csv" --resume
```

## Building Models

### Required Inputs (Model Building)

`--assemblypath` : Path to your assembly directory.

- **Must** be a full path, not a relative one.
- Assemblies **must** have .fasta extension (otherwise the pipeline breaks)

`--hostdata` : Path to your metadata file. The following fields must be present:

- **Filename**: Full filename of your assembly, including the extension (.fasta)
- **Host**: Host your assembly was obtained from. Can be as many as you want or as few as two. This pipeline assumes you have human hosts in your dataset. If not, the R-script modules which generate human-scoring models may crash.
- **Region**: Country/Region the assembly was obtained from. It is used for the detection for clonal clusters. Fields can be left blank, but must be present in the metadata file.
- **Year:** Year the assembly was collected. It is used for the detection for clonal clusters. Fields can be left blank, but must be present.
- Below is an example of the contents of the metadata file:
    ```
    Filename,Source.Host,Year,Region
    SAL_AB4979AA_AS.result.fasta,Bovine,2001,USA
    SAL_AB8755AA_AS.result.fasta,Bovine,2001,USA
    SAL_BA6077AA_AS.scaffold.fasta,Bovine,2001,USA
    ```
- An alternative minimal metadata file with the Year and Region fields blank:
    ```
    Filename,Source.Host,Year,Region
    SAL_AB4979AA_AS.result.fasta,Bovine,,
    SAL_AB8755AA_AS.result.fasta,Bovine,,
    SAL_BA6077AA_AS.scaffold.fasta,Bovine,,
    ```

Please refer to the metadata.csv file inside the input test folder for another example.

### Advanced Runs (Model Building)

- Specify output directory to ./pipeline_results instead of the default ./out :

``` bash
nextflow run pipeline_dsl2.nf --assemblypath "/home/username/myproject/input_data" --hostdata "/home/username/myproject/input_data/metadata.csv" --outdir "./pipeline_results"
```

- Metadata file with different headings (eg a metdata file that has the headers File,Source,Area,Collection.Year instead of the expected Filename,Source.Host,Region,Year):

``` bash
nextflow run pipeline_dsl2.nf --assemblypath "/home/username/myproject/input_data" --hostdata "/home/username/myproject/input_data/metadata.csv" --assembly_column="File" --host_column="Host" --region_columns="Area" --year_collection="Collection.Year"
```

- Change the cutoffs for filtering out low quality assemblies:

``` bash
nextflow run pipeline_dsl2.nf --assemblypath "/home/username/myproject/input_data" --hostdata "/home/username/myproject/input_data/metadata.csv" --as_ln_upr=7000000 --as_ln_lwr=3000000 --ctg_count=100 --largest_ctg=150000 --n50=30000 --gc_upr=56 --gc_lwr=40
```

- Change the threshold for SNP difference used to cluster assemblies for clonal detection:

``` bash
nextflow run pipeline_dsl2.nf --assemblypath "/home/username/myproject/input_data" --hostdata "/home/username/myproject/input_data/metadata.csv" --snp_dist_threshold=100
```

- Provide your own trusted protein file for prokka-based annotations:

``` bash
nextflow run pipeline_dsl2.nf --assemblypath "/home/username/myproject/input_data" --hostdata "/home/username/myproject/input_data/metadata.csv" --prokka_ref "./my_folder/protein_file.fasta"
```

- Change amrfinder species & have prokka use the Escherichia genus for annotation:

``` bash
nextflow run pipeline_dsl2.nf --assemblypath "/home/username/myproject/input_data" --hostdata "/home/username/myproject/input_data/metadata.csv" --amr_species="ecoli" --prokka_extra "--genus Escherichia"
```

- Run panaroo on a different filtering setting:

``` bash
nextflow run pipeline_dsl2.nf --assemblypath "/home/username/myproject/input_data" --hostdata "/home/username/myproject/input_data/metadata.csv"--panaroo_mode="strict"
```

- Provide a SNP reference file:

``` bash
nextflow run pipeline_dsl2.nf --assemblypath "/home/username/myproject/input_data" --hostdata "/home/username/myproject/input_data/metadata.csv"  --snp_ref "./my_folder/my_snp_ref.fasta"
```

- Change the number of cores for panaroo and piggy to run on:

``` bash
nextflow run pipeline_dsl2.nf --assemblypath "/home/username/myproject/input_data" --hostdata "/home/username/myproject/input_data/metadata.csv"  --threads=5
```

TODO Run with docker

TODO Run with singularity

TODO add -resume and profile in paarmeters?

### Full Parameters (Model Building)

To get a full list of the available parameters, run: `nextflow  run pipeline_dsl2.nf --help`

``` text
Usage:
        The options for running the pipeline to build models are structured as follows:
        --option                       Description/Notes [default value]

    Mandatory Parameters:
         --assemblypath                Directory of your fasta files (full path required) [./test_data/model_build_input_test]
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
        --as_ln_lwr                    Minimum accepted assembly length [4000000]
        --ctg_count                    Minimum accepted number of contigs [500]
        --largest_ctg                  Minimum accepted length of largest contig [100000] 
        --n50                          Minimum accepted n50 [50000]
        --gc_upr                       Minimum accepted GC % [54]
        --gc_lwr                       Maximum accepted GC % [50]

    Clonal Detection Parameters:
        --snp_dist_threshold           SNP difference used to detect clonal clusters. Used in conjunction with region & collection year metadata. [10]
        
    Tool-specific Parameters:
        --prokka_ref                   'Trusted' protein file for prokka (prokka --proteins) [./data/stm_proteinref.fasta]
        --amr_species                  Assembly species for amrfinder (amrfinder -O) ["Salmonella"]
        --panaroo_mode                 Panaroo assembly filtering mode (panaroo --clean-mode) ["moderate"]
        --snp_ref                      Reference file for snippy (snippy --ref) [./data/stm_sl1344.fasta]
        --threads                      Num of threads to use for panaroo & piggy (panaroo -t & piggy -t) [10]
    
    If you wish to alter the scripts used to generate the models, simply edit the appropriate 'model_building' R scripts in ./data/ - ONLY DO SO IF YOU KNOW WHAT YOU ARE DOING
```

### Outputs (Model Building)

The output folder should contain the following folders:

- `0.Reports/`: Contains housekeeping  files about the pipeline execution, which include:
  - [Execution report](https://www.nextflow.io/docs/latest/tracing.html#execution-report) as `report.html`
  - [Trace Report](https://www.nextflow.io/docs/latest/tracing.html#trace-report) as `trace.txt`
  - [DAG Visualization](https://www.nextflow.io/docs/latest/tracing.html#dag-visualisation) as `DAG.svg`
  - [Timeline Report](https://www.nextflow.io/docs/latest/tracing.html#timeline-report) as `timeline.html`
  - `1.assembly_quality`
    - `good_noext_metadata.csv`
    - `good_assemblies`
      - `*.fasta`
- `2.genomic_features/`
  - `amr_all.tsv`
  - `scoary_traitfile.csv`
  - `annotations/`
    - `*.gff`
    - `*.faa`
    - `*.fna`
  - `pv_out`
    - `panaroo_out/`
      - `gene_presence_absence.Rtab`
      - `pan_genome_reference.fa`
    - `roary_out/`
      - `gene_presence_absence_roary.csv`
    - `scoary_pv/`
      - `*.results.csv`
  - `igr_out/`
    - `piggy_out/`
      - `IGR_presence_absence.csv`
      - `IGR_presence_absence.Rtab`
      - `representative_clusters_merged.fasta`
    - `scoary_igr/`
      - `*.results.csv`
  - `snp_core_out/`
    - `core.tab`
    - `core.aln`
  - `3.clonal_detection/`
    - `base_cluster_static.png`
    - `good_nonclonal_metadata.csv`
    - `clusters/`
      - `*.list`
- `4.model/`
  - `model_ref/`
    - `igr_filter.fasta`
    - `pv_filter.fasta`
    - `assemblyblastdb*`
  - `model_input/`
    - `amr_gene_all.tsv`
    - `amr_gene_bps.tsv`
    - `amr_gene_human.tsv`
    - `amr_class_all.tsv`
    - `amr_class_bps.tsv`
    - `amr_class_human.tsv`
    - `igr_all.tsv`
    - `igr_bps.tsv`
    - `igr_human.tsv`
    - `pv_all.tsv`
    - `pv_bps.tsv`
    - `pv_human.tsv`
    - `snp_abudance_all.tsv`
    - `snp_abudance_bps.tsv`
    - `snp_abudance_human.tsv`
  - `models/`
    - `*.rds`
  - `predictions/`
    - `*.csv`
  - `plots/`
    - `*.png`

TODO ADD Descriptions of output files

### Workflow (Model Building)

![Simplified Diagram of Model Building Pipeline](https://github.com/Antonia-Chalka/stm_ml_pipeline/blob/main/data/pipeline_diagrams-build.png?raw=true)

### Clonal Filtering

In order to maximize genetic diversity and reduced bias resulting from overrepresented features as a result of clonal outbreaks, the pipeline has a clonal filtering stage. Assemblies are placed in a cluster if **all** the following conditions are met:

- Same host
- Same region
- Same year
- <10 SNP difference (threshold can be changed by the `--snp_dist_threshold` flag)

Assemblies that fullfil all of the above conditions are collected into clonal clusters. **One** assembly is chosen in random form each cluster and the rest are removed from the training dataset.

## Testing Models

### Required Inputs (Testing Models)

TODO WRITEUP

### Advanced Runs (Testing Models)

TODO Include in write up that metadata folder can just be filenames

TODO WRITEUP

### Full Parameters (Testing Models)

TODO WRITEUP

### Outputs (Testing Models)

TODO WRITEUP

### Workflow (Testing Models)

![Simplified Diagram of Model Testing Pipeline](https://github.com/Antonia-Chalka/stm_ml_pipeline/blob/main/data/pipeline_diagrams-test.png?raw=true)

## Benchmarks

**Specs**:

- CPU: AMD Ryzen 9 5900X (12 cores, 24 logical processors)

- 128 GB RAM

**Test Dataset:**

``` 
Completed at: 09-Dec-2021 19:00:52
Duration    : 17m 6s
CPU hours   : 3.7
Succeeded   : 80
```

TODO update when u fix test dataset issue

**Large Dataset (~5k STm sequences):**

TODO add benchmarks for big dataset

## Citation

TODO ADD CITATION
