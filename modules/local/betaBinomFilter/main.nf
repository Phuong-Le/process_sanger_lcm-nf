process betaBinomFilter {
    publishDir "${params.outdir}/filter_${mut_type}_out/${pdid}", overwrite: true
    
    // filtering according to the indices obtained from betaBinomFilterIndex
    input:
    tuple val(pdid), val(sample_id), val(match_normal_id), path(vcf_to_filter), path(vcf_tbi_to_filter), path(bed_idx)
    val mut_type

    output:
    path filtered_vcf

    script:
    vcf_names = vcf_to_filter.getName().tokenize(".")
    vcf_extension = vcf_names.tail().init().join(".") // removing the filename and the gz extension
    filtered_vcf="${sample_id}.somatics_filtered.artefacts_filtered.${vcf_extension}"
    """
    tabix -h -R ${bed_idx} ${vcf_to_filter} > ${filtered_vcf}
    """
}