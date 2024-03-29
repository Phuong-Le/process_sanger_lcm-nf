process pindelFilter {
    publishDir "${params.outdir}/${pdid}", overwrite: false, pattern: "*.filter.vcf.gz"
    
    // hairpin filter the mutation file 
    input:
    tuple val(sample_id), val(match_normal_id), val(pdid), path(vcf), path(bam), path(bai), path(bam_match), path(bai_match)
    path vcfilter_config

    output: 
    tuple val(pdid), val(sample_id), val(match_normal_id), path("*.filter.vcf.gz"), path("*.filter.vcf.gz.tbi"), path(bam), path(bai), path(bam_match), path(bai_match)

    script:
    """
    pindel_filter.sh --vcf ${vcf} --vcfilter_config ${vcfilter_config} 
    """
}