process concatMatrices {


    input:
    path matrix_dirs, stageAs: "indir/*"

    output:
    path outdir

    script:
    outdir = "combined_matrices_by_branches"
    """
    mkdir ${outdir}
    concat_mutmats.py --indir indir --outdir ${outdir}
    """
}