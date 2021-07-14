#!/usr/bin/env nextflow

params.assemblypath = "$baseDir/input_test/"
params.hostdata = "$baseDir/input_test/metadata.csv"
params.outdir = "$baseDir/out/"

assemblies = Channel.fromPath("${params.assemblypath}/*.fasta")

params.as_ln_upr = 6000000
params.as_ln_lwr = 4000000
params.ctg_count = 500
params.largest_ctg = 100000
params.n50 = 50000
params.gc_upr = 54
params.gc_lwr = 50

process assembly_qc {
    
    input:
    file assembly from assemblies

    output:
    file "${assembly.baseName}_qc/transposed_report.tsv" into quast_ch

    script:
    """
    quast.py  -o "${assembly.baseName}_qc" "${assembly}" 
    """
}

process printqc {
    publishDir  "${params.outdir}/qc_report", mode: 'copy', overwrite: true

    input:
    file qcs from quast_ch.collectFile(keepHeader:true, skip:1)

    output:
    file 'pass_qc.tsv' into good_qc

    script:
    """
    awk -F "\t" '{ if((\$2 < $params.ctg_count) && (\$8 > $params.as_ln_lwr && \$8 < $params.as_ln_upr) && (\$15 > $params.largest_ctg) && (\$17 > $params.gc_lwr && \$17 < $params.gc_upr) && (\$18 > $params.n50)) { print } }' $qcs > pass_qc.tsv

    mkdir ${params.outdir}good_assemblies/ -p
    
    for assembly_name in `cut -f1 pass_qc.tsv`
    do 
        echo \$assembly_name
        cp ${params.assemblypath}/\${assembly_name}.fasta ${params.outdir}/good_assemblies/\${assembly_name}.fasta
    done
    """

}

good_assemblies_ch = Channel.fromPath("${params.outdir}/good_assemblies/*.fasta")

