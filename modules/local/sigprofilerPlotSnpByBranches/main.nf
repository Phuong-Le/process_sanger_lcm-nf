process sigprofilerPlotSnpByBranches {
    publishDir "${params.outdir}/${outdir_basename}", overwrite: false

    input:
    path mutmat_dir
    val outdir_basename

    output:
    path mutmat_plot_dir

    script:
    mutmat_plot_dir = "combined_matrices_by_branches/plots"
    """
    sigprofiler_plotting.py --matrix_dir ${mutmat_dir} --project branch_combined_mutmat --output_path ${mutmat_plot_dir}
    """
}