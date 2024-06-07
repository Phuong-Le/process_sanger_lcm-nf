process sigprofilerPlotSnp {
    publishDir "${params.outdir}", overwrite: false

    input:
    path mutmat_dir, stageAs: 'sample_mutmat/*'

    output:
    path mutmat_plot_dir

    script:
    mutmat_plot_dir = "sample_mutmat/plots"
    """
    sigprofiler_plotting_snp.py --matrix_dir ${mutmat_dir} --project sample_mutmat --output_path ${mutmat_plot_dir}
    """
}