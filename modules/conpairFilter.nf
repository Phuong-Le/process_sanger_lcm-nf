process conpairFilter {
    publishDir "${params.outdir}/conpair_out", overwrite: false

    input:
    path concordance
    path contamination
    path samples_path


    output:
    path outfile
    path "*.log"
    path concordance
    path contamination

    script:
    outfile = "sample_paths_contamination_filtered.txt"
    """
    conpair_contamination_filter.py --samples_path ${samples_path} --concordance_path ${concordance} --contamination_path ${contamination} --concordance_threshold ${params.concordance_threshold} --contamination_threshold_samples ${params.contamination_threshold_samples} --contamination_threshold_match ${params.contamination_threshold_match} --outfile ${outfile}
    """
}