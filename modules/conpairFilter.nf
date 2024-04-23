process conpairFilter {
    publishDir "${params.outdir}", overwrite: false

    input:
    path concordance
    path contamination
    path sample_paths


    output:
    path outfile
    path "conpair_filter.log"

    script:
    outfile = "sample_paths_new.txt"
    """
    conpair_filter.py --concordance $concordance --contamination $contamination --sample_paths $sample_paths --outfile $outfile 
    """
}