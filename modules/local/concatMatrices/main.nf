process concatMatrices {
    publishDir "${params.outdir}/${outdir_basename}", overwrite: false

    input:
    path matrix_dirs, stageAs: "indir/matrix_generator*"
    val outdir_basename

    output:
    path outdir

    script:
    outdir = "combined_matrices_by_branches"
    """
    mkdir ${outdir}
    concat_mutmats.py --indir indir --outdir ${outdir}
    """
}