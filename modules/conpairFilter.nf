process conpairFilter {
    publishDir "${params.outdir}", overwrite: false

    input:
    path concordance
    path contamination
    path samples_path


    output:
    path outfile
    path "*.log"

    script:
    outfile = "sample_paths_contamination_filtered.txt"
    """
    conpair_contamination_filter.py --samples_path ${samples_path} --concordance_path ${concordance} --contamination_path ${contamination} --concordance_threshold ${params.concordance_threshold} --contamination_threshold ${params.contamination_threshold} --outfile ${outfile}
    """
}