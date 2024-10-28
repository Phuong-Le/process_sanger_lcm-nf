process matrixGeneratorOnSamples {
    publishDir "${params.outdir}/filter_${mut_type}_out", overwrite: false

    input:
    path vcf_ls, stageAs: 'sample_mutmat/*'
    val mut_type

    output:
    path mutmat_dir

    script:
    mutmat_dir = "sample_mutmat/output"
    """
    SigProfilerMatrixGenerator matrix_generator sample_mutmat GRCh38 sample_mutmat
    """
}