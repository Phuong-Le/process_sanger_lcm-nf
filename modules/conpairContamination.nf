process conpairContamination {
    input:
    tuple val(match_normal_id), path(match_pileup), val(sample_id), path(sample_pileup)

    output:
    path(contamination_path)

    script:
    contamination_path = "${sample_id}_${match_normal_id}.contamination.txt"
    """
    estimate_tumor_normal_contamination.py -T ${sample_pileup} -N ${match_pileup} -M ${params.marker_txt} -O ${contamination_path}
    """

}