process getPhylogeny {
    publishDir "${params.outdir}/${outdir_basename}/${pdid}", overwrite: true

    input:
    tuple val(pdid), path(nr_path), path(nv_path), path(genotype_bin_path)
    val outdir_basename

    output:
    tuple val(pdid), path(phylogeny_vcf)
    tuple val(pdid), path(tree_path)
    tuple val(pdid), path(fasta), path(tree_w_branchlength), path(tree_w_branchlength_plot)
    path mpboot_log

    script:
    fasta = "${pdid}.fasta"
    tree_path = "${fasta}.treefile"
    mpboot_log = "${fasta}.log"
    tree_w_branchlength = "${pdid}.tree_with_branch_length.tree"
    tree_w_branchlength_plot = "${pdid}.tree_with_branch_length.pdf"
    phylogeny_vcf = "${pdid}.muts_assigned_to_tree.txt"
    """
    genotype_bin_to_fasta.py --genotype_bin_path ${genotype_bin_path} --outfile ${fasta}
    mpboot -s ${fasta} -bb 1000 -cost e -st DNA
    Rscript --vanilla ${projectDir}/bin/assign_mut_to_tree.R --nr_path=$nr_path --nv_path=$nv_path --genotype_bin_path=$genotype_bin_path --tree_path=$tree_path --out_prefix=$pdid
    """

}