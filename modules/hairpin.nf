process hairpinFilter {
    publishDir "${params.outdir}/${pdid}", overwrite: false

    // hairpin filter the mutation file 
    input:
    tuple val(sample_id), val(match_normal_id), val(pdid), path(vcf), path(vcf_tbi), path(bam), path(bai), path(bas), path(met), path(bam_match), path(bai_match)
    path vcfilter_config

    output: 
    tuple val(pdid), val(sample_id), val(match_normal_id), path("*.hairpin.filter.vcf.gz"), path("*.hairpin.filter.vcf.gz.tbi"), path(bam), path(bai), path(bam_match), path(bai_match)

    script:
    """
    hairpin.sh --vcf ${vcf} --bam ${bam} --vcfilter_config ${vcfilter_config} 
    """
}