process matrixGeneratorOnSamples {
    publishDir "${params.outdir}", overwrite: false

    input:
    path vcf_ls, stageAs: 'sample_mutmat/*'

    output:
    path publish_dir

    script:
    publish_dir = "sample_mutmat/output"
    """
    SigProfilerMatrixGenerator matrix_generator sample_mutmat GRCh38 sample_mutmat
    """
}