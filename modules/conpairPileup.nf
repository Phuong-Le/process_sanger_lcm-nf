process conpairPileup {
    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(pileup)

    script:
    """
    """

}