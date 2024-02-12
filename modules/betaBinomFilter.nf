process betaBinomFilter {
    publishDir "${params.outdir}/${pdid}", overwrite: false
    
    // filtering according to the indices obtained from betaBinomFilterIndex
    input:
    tuple val(pdid), val(sample_id), val(match_normal_id), path(vcf_to_filter), path(vcf_tbi_to_filter), path(bed_idx)

    output:
    tuple val(sample_id), path(filtered_vcf)
    
    script:
    vcf_names = vcf_to_filter.getName().tokenize(".")
    vcf_extension = vcf_names.tail().join(".")
    filtered_vcf="${sample_id}.${bed_idx.getSimpleName()}.${vcf_extension}"
    """
    tabix -h -R ${bed_idx} ${vcf_to_filter} > ${filtered_vcf}
    """
}