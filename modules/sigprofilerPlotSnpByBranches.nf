process sigprofilerPlotSnpByBranches {
    publishDir "${params.outdir}", overwrite: false

    input:
    path mutmat_dir

    output:
    path mutmat_plot_dir

    script:
    mutmat_plot_dir = "combined_matrices_by_branches/plots"
    """
    sigprofiler_plotting_snp.py --matrix_dir ${mutmat_dir} --project branch_combined_mutmat --output_path ${mutmat_plot_dir}
    """
}