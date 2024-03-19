process conpairPileup {
    input:
    tuple val(match_normal_id), val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path(pileup)

    script:
    pileup = "${sample_id}.pileup"
    """
    gatk --java-options -Xmx10g Pileup -R ${params.reference_genome} -I ${bam} -L ${params.marker_bed} -O ${pileup} -verbose -RF NotDuplicateReadFilter -RF CigarContainsNoNOperator -RF MatchingBasesAndQualsReadFilter
    """

}