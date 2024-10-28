process mutToTree {
    publishDir "${params.outdir}/${outdir_basename}/${pdid}", overwrite: true

    input:
    tuple val(pdid), path(topology), path(nr_path), path(nv_path), path(genotype_bin_path)
    val outdir_basename

    output:
    tuple val(pdid), path(phylogeny_vcf_with_header)
    tuple val(pdid), path(tree_w_branchlength), path(tree_w_branchlength_plot)

    script:
    tree_w_branchlength = "${pdid}.tree_with_branch_length.tree"
    tree_w_branchlength_plot = "${pdid}.tree_with_branch_length.pdf"
    phylogeny_vcf_with_header = "${pdid}.muts_assigned_to_tree.txt"
    """
    Rscript --vanilla ${projectDir}/bin/assign_mut_to_tree.R --nr_path=$nr_path --nv_path=$nv_path --genotype_bin_path=$genotype_bin_path --tree_path=$topology --out_prefix=$pdid
    """

}