process conpairPileup {
    publishDir "${params.outdir}/conpair_out/pileup", overwrite: false

    input:
    tuple val(match_normal_id), val(sample_id), path(bam), path(bai)
    path(marker_bed)
    path(reference_genome)
    path(reference_genome_dict)
    path(reference_genome_idx)

    output:
    tuple val(match_normal_id), val(sample_id), path(pileup)

    script:
    pileup = "${sample_id}.pileup"
    """
    gatk --java-options -Xmx10g Pileup -R ${reference_genome} -I ${bam} -L ${marker_bed} -O ${pileup} -verbose -RF NotDuplicateReadFilter -RF CigarContainsNoNOperator -RF MatchingBasesAndQualsReadFilter
    """

}