process betaBinomFilterIndex {
    publishDir "${params.outdir}/${params.mut_type}/${pdid}", overwrite: false

    // Beta Binomial filtering of germline mutations and artefacts, based on Tim Coorens' R script
    // The outcome is a bed file for the indices of the PASSED mutations
    input:
    tuple val(pdid), val(sample_id_ls), val(match_normal_id), path(vaf)

    output:
    tuple val(pdid), path("*.bed")
    path("germline_ids.txt") 
    path("somatic_ids.txt")
    path("somatic_ids_rho.txt")


    script:
    """
    Rscript --vanilla ${projectDir}/bin/beta_binom_filter_index.R --libpath=${projectDir}/lib --cgpvaf_out=${vaf} --match_normal_id=${match_normal_id} --outdir=. 
    """

}