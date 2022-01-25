# Workflow for Host Attribution Machine Learning Models 

A DSL2 Nextflow & Docker pipeline used to build bacterial source attribution machine learning models from assembled genomes, using SNPs, protein variants (pvs), intergenic regions (igrs) and AMR profiles. Currently, it has only been tested for *Salmonella typhimurium* sequences.

As used in:
TODO Add citation

## Installation

* [Install Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)

* [Install Docker](https://docs.docker.com/get-docker/)

* Pull Required Docker/Singularity Images

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
docker pull annitachalka/r_model_build:1.01

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
```

## Usage

Simplest way to run:

`nextflow run pipeline_dsl2.nf --assemblypath "/home/username/myproject/input_data" --hostdata "/home/username/myproject/input_data/metadata.csv" `

More advanced run:

`nextflow run pipeline_dsl2.nf --assemblypath "/home/username/myproject/input_data" --hostdata "/home/username/myproject/input_data/metadata.csv" --snp_dist_threshold=100 --year_collected="Year.Obtained" --panaroo_mode='strict` -resume

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

## Additional Parameters/Inputs
 TODO ADD OPTIONAL PARAMETER GUIDE
  https://github.com/AdmiralenOla/Scoary/blob/master/README.md

`--outdir`

`--assembly_column`
`--host_column`
`--region_column`
`--year_collected`

`--as_ln_upr` = 6000000
`--as_ln_lwr` = 4000000
`--ctg_count` = 500
`--largest_ctg` = 100000
`--n50` = 50000
`--gc_upr` = 54
`--gc_lwr` = 50

`--snp_dist_threshold`

`--prokka_ref` = "$projectDir/data/stm_proteinref.fasta" 
`--snp_ref` = "$projectDir/data/stm_sl1344.fasta"

`--amr_species`='Salmonella'
`--panaroo_mode` = 'moderate'

## Outputs


## Workflow Outline

TODO ADD DIAGRAM

## Benchmarks

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
