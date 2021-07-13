#!/usr/bin/env nextflow

params.assemblypath = "$baseDir/input_test/*.fasta"
params.hostdata = "$baseDir/input_test/metadata.csv"

assemblies = Channel.fromPath(params.assemblypath)

process assembly_qc {
    
    input:
    file assembly from assemblies

    output:
    file "${assembly.baseName}_qc/report.tsv" into quast_ch

    """
    quast.py  -o "${assembly.baseName}_qc" "${assembly}" 
    """
}

