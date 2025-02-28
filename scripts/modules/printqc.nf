// Filter assemblies based on quast-derived metrics
process printqc {
    publishDir  "${params.outdir}/1.assembly_quality/good_assemblies/", mode: 'copy', overwrite: true, pattern : "*.${params.fileextension}"
    publishDir  "${params.outdir}/1.assembly_quality/", mode: 'copy', overwrite: true, pattern : 'good_noext_metadata.csv'

    cache 'lenient'

    input:
    path qcs
    path metadata_file

    output:
    path "*.${params.fileextension}", emit: good_assemblies
    path 'good_noext_metadata.csv', emit: good_metadata
    path 'pass_qc.tsv', emit: good_assemblies_list

    script:
    """
    # Set error handling
    set -e
    
    # Check if input files & assembly path exist
    if [ ! -f "${qcs}" ]; then
        echo "ERROR: QC file "${qcs}" not found"
        exit 1
    fi

    if [ ! -f "${metadata_file}" ]; then
        echo "ERROR: Metadata file "${metadata_file}" not found"
        exit 1
    fi

    if [ ! -d "${params.assemblypath}" ]; then
        echo "ERROR: Assembly directory "${params.assemblypath}" not found"
        exit 1
    fi

    # Create list of assemblies that pass qc & Copy them
    awk -F "\t" '{ 

        if( (\$2 < ${params.ctg_count}) && 
            (\$8 > ${params.as_ln_lwr} && \$8 < ${params.as_ln_upr}) && 
            (\$15 > ${params.largest}_ctg) && 
            (\$17 > ${params.gc_lwr} && \$17 < ${params.gc_upr}) && 
            (\$18 > ${params.n50})) { 
            print 
        } 
    }' "${qcs}" > pass_qc.tsv

    # Check if any assemblies passed QC
    if [ ! -s pass_qc.tsv ]; then
        echo "ERROR: No assemblies passed QC criteria"
        exit 100
    fi

    # Create symlinks for passing assemblies
    while IFS= read -r line; do
        assembly_name=\$(echo "\$line" | cut -f1)
        assembly_path="${params.assemblypath}/\${assembly_name}.${params.fileextension}"
        
        # Check if source file exists before creating symlink
        if [ -f "\$assembly_path" ]; then
            ln -sf "\$assembly_path" "\${assembly_name}.${params.fileextension}"
        else
            echo "WARNING: Assembly file '\$assembly_path' not found, skipping"
        fi
    done < pass_qc.tsv
    
    # Add headers to new metadata file by copying the first line of the orginal metadata file
    head -n 1 "${metadata_file}" > good_metadata.csv
    
    # Create list of good assemblies that were successfully linked
    ls -1 *.${params.fileextension} 2>/dev/null > good_assemblies.txt || touch good_assemblies.txt
    
    # Append matching metadata entries
    while IFS= read -r assembly_file; do
        assembly_name=\$(basename "\$assembly_file" .${params.fileextension})
        grep -m 1 -F "\$assembly_name" "${metadata_file}" >> good_metadata.csv || echo "WARNING: No metadata found for \$assembly_name"
    done < good_assemblies.txt


    # Create version with extensions removed
    awk 'NR==1 {print; next} 
        { gsub(/\\.${params.fileextension}/,"", \$1); print }' good_metadata.csv > good_noext_metadata.csv
    
    """
}
