#!/usr/bin/env nextflow

params.assemblypath = "$baseDir/input_test/*.fasta"
params.hostdata = "$baseDir/input_test/metadata.csv"

assemblies = Channel.fromPath(params.assemblypath)

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

    """
    quast.py  -o "${assembly.baseName}_qc" "${assembly}" 
    """
}

process printqc {
    input:
    file qcs from quast_ch.collectFile(keepHeader:true, skip:1)

    output:
    file 'pass_qc.tsv' into good_qc
    script:
    """
    awk -F "\t" '{ if((\$2 < $params.ctg_count) && (\$8 > $params.as_ln_lwr && \$8 < $params.as_ln_upr) && (\$15 > $params.largest_ctg) && (\$17 > $params.gc_lwr && \$17 < $params.gc_upr) && (\$18 > $params.n50)) { print } }' $qcs > pass_qc.tsv

    cd $params.assemblypath
    ls
    """

}

good_qc
    .splitCsv(header: false, sep: "\t" )
    .map{ row-> row[0] }
    .view()


