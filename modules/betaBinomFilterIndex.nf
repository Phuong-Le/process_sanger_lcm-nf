process betaBinomFilterIndex {
    publishDir "${params.outdir}/${pdid}", overwrite: false

    // Beta Binomial filtering of germline mutations and artefacts, based on Tim Coorens' R script
    // The outcome is a bed file for the indices of the PASSED mutations
    input:
    tuple val(pdid), val(sample_id_ls), val(match_normal_id), path(vaf)

    output:
    tuple val(pdid), path("*.bed")
    path "germline_ids.txt", optional: true 
    path "somatic_ids.txt", optional: true 
    path "somatic_ids_rho.txt", optional: true 
    tuple val(pdid), path("NR_somatic_noartefacts.txt"), path("NV_somatic_noartefacts.txt"), path("genotype_bin.txt"), optional: true

    script:
    """
    Rscript --vanilla ${projectDir}/bin/beta_binom_filter_index.R --libpath=${projectDir}/lib --cgpvaf_out=${vaf} --match_normal_id=${match_normal_id} --outdir=. 
    """

}