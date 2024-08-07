process matrixGeneratorOnBranches {
    publishDir "${params.outdir}/${outdir_basename}/${pdid}", overwrite: false

    input:
    tuple val(pdid), path(phylogeny_vcf_with_header)
    val outdir_basename

    output:
    path matrix_generator_dir

    script:
    outdir = 'matrix_by_branch'
    matrix_generator_dir = "${outdir}/matrix_generator"
    """
    mkdir ${outdir}
    split_vcf_to_branch.py --vcf_path ${phylogeny_vcf_with_header} --outdir ${outdir} --prefix ${pdid}
    SigProfilerMatrixGenerator matrix_generator ${pdid} GRCh38 ${matrix_generator_dir}
    rm ${matrix_generator_dir}/*.vcf
    """
}