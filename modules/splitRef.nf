process splitRef {
    publishDir "${params.reference_genome_cachedir}", mode: "copy", overwrite: false
    
    // rename the reference files by its identifier 
    input:
    path fasta

    output:
    path "*.fa"

    script:
    """
    split_ref.sh --reference ${fasta}
    """

}