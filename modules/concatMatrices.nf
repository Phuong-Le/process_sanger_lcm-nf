process concatMatrices {
    publishDir "${params.outdir}", overwrite: false

    input:
    path matrix_dirs, stageAs: "indir/matrix_generator*"

    output:
    path outdir

    script:
    outdir = "combined_matrices_by_branches"
    """
    mkdir ${outdir}
    concat_mutmats.py --indir indir --outdir ${outdir}
    """
}