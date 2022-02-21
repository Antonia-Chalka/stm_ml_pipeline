# Workflow for Host Attribution Machine Learning Models 

A DSL2 Nextflow & Docker pipeline used to build bacterial source attribution machine learning models from assembled genomes, using SNPs, protein variants (PVs), intergenic regions (IGRs) and AMR profiles. Currently, it has only been tested for *Salmonella typhimurium* sequences.

As used in:

TODO Add citation

# Installation

* [Install Java](https://java.com/en/download/help/download_options.html)

* [Install Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)

* [Install Docker](https://docs.docker.com/get-docker/)

* Pull Required Docker/Singularity Images

  * With Docker:

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

  * With SIngularity:
  
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

# Building Models

## Usage

Simplest way to run:

`nextflow run pipeline_dsl2.nf --assemblypath "/home/username/myproject/input_data" --hostdata "/home/username/myproject/input_data/metadata.csv"`

More advanced run:

`nextflow run pipeline_dsl2.nf --assemblypath "/home/username/myproject/input_data" --hostdata "/home/username/myproject/input_data/metadata.csv" --snp_dist_threshold=100 --year_collected="Year.Obtained" --panaroo_mode='strict -resume`

TODO Add more advanced runs

* Metadata file
* Assembly stats
* Threads
* Organism
* resume option

TODO Add testing runs

* From previous run
* From importing models

## Required Inputs

`--assemblypath` : Folder where your assemblies are kept.

* **Must** be a full path, not a relative one.

* Assemblies **must** have .fasta extension (otherwise the pipeline breaks)

`--hostdata` : File where your metadata is places. The following fields must be present:

* **Filename**: Full filename of your assembly, along with the extsnion (.fasta)

* **Host**: Host your assembly was obtained from. Can be as many or as few as two. This pipeline assumes you have human hosts in your dataset. If not, the later R-script modules may crash.

* **Region**: Country/Region the assembly was obtained from. It is used for the detection for clonal clusters. Fields can be left blank, but must be present.

* **Year:** Year the assembly was collected. It is used for the detection for clonal clusters. Fields can be left blank, but must be present.

Please refer to the metadata.csv file inside the input test folder for an example.

## Full Parameters/Inputs

To get a full list of the available paramteters, run: `nextflow  run pipeline_dsl2.nf --help`

``` text
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
```
## Outputs
The output folder should contain the following folders: 

0.Reports: Contains housekeeping  files about the pipeline execution, which include:

* [Execution report](https://www.nextflow.io/docs/latest/tracing.html#execution-report) as report.html
* [Trace Report](https://www.nextflow.io/docs/latest/tracing.html#trace-report) as trace.txt
* [DAG Visualisation](https://www.nextflow.io/docs/latest/tracing.html#dag-visualisation) as DAG.svg
* [Timeline Report](https://www.nextflow.io/docs/latest/tracing.html#timeline-report) as timeline.html

## Workflow

![Simplified Diagram of Model Building Pipeline](https://github.com/Antonia-Chalka/stm_ml_pipeline/blob/main/data/STm_plan.png?raw=true)

## Assembly QC

TODO WRITEUP

## Clonal Filtering

TODO WRITEUP


TODO ADD DIAGRAM

TODO WRITEUP

# Testing Models

TODO WRITEUP

## Usage

TODO WRITEUP

## Required Input

## Full Parameters

TODO WRITEUP

## Output

TODO WRITEUP

# Workflow

## Building Models Pipeline

TODO ADD DAG FILE

## Testing Models Pipeline

TODO ADD FILE

# Benchmarks

**Specs**:

* CPU: AMD Ryzen 9 5900X (12 cores, 24 logical processors)

* 64 GB RAM

Input test data (~12 final dataset assemblies):

``` 
Completed at: 09-Dec-2021 19:00:52
Duration    : 17m 6s
CPU hours   : 3.7
Succeeded   : 80
```

My STm Dataset (3313 assemblies that pass qc, ~2.2k nonclonal):
TODO add benchmarks for big dataset
