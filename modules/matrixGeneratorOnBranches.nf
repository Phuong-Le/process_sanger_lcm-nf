process matrixGeneratorOnBranches {
    publishDir "${params.outdir}/${pdid}", overwrite: false

    input:
    tuple val(pdid), path(phylogeny_vcf)

    output:
    path matrix_generator_dir
    path vcf_with_header

    script:
    outdir = 'matrix_by_branch'
    matrix_generator_dir = "${outdir}/matrix_generator"
    vcf_with_header = "${outdir}/vcf_with_header"
    """
    mkdir ${outdir}
    split_vcf_to_branch.py --vcf_path ${phylogeny_vcf} --outdir ${outdir} --prefix ${pdid}
    SigProfilerMatrixGenerator matrix_generator ${pdid} GRCh38 ${matrix_generator_dir}
    rm ${matrix_generator_dir}/*.vcf
    """
}