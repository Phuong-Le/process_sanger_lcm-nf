process betaBinomFilter {
    publishDir "${params.outdir}/${pdid}", overwrite: false
    
    // filtering according to the indices obtained from betaBinomFilterIndex
    input:
    tuple val(pdid), val(sample_id), val(match_normal_id), path(vcf_to_filter), path(vcf_tbi_to_filter), path(bed_idx)

    output:
    tuple val(sample_id), path(final_vcf)
    path filtered_vcf

    script:
    vcf_names = vcf_to_filter.getName().tokenize(".")
    vcf_extension = vcf_names.tail().init().join(".") // removing the filename and the gz extension
    filtered_vcf="${sample_id}.somatics_filtered.artefacts_filtered.${vcf_extension}"
    final_vcf="${sample_id}.final.vcf"
    """
    tabix -h -R ${bed_idx} ${vcf_to_filter} > ${filtered_vcf}
    grep -vE '^(##)' ${filtered_vcf} | sed 's/^#//' > ${final_vcf}
    """
}