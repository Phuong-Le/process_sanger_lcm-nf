process getMutMat {

    // get mutation matrix from vcf
    input:
    tuple val(sample_id), path(final_vcf)
    path reference
    val kmer

    output:
    path mutmat_path

    script:
    mutmat_path = "${sample_id}.txt"
    """
    get_mutation_matrix.py --ref_dir . --vcf_path ${final_vcf} --sample_id ${sample_id} --outfile ${mutmat_path} --kmer ${kmer}
    """
}   