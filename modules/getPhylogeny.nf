process getPhylogeny {
    input:
    tuple val(pdid), path(nr_path), path(nv_path), path(genotype_bin_path)

    output:
    tuple val(pdid), path(fasta), path(tree_path), path(tree_w_branchlength), path(tree_w_branchlenth_plot), path(phylogeny_vcf)
    path mpboot_log

    script:
    fasta = "${pdid}.fasta"
    tree_path = "${fasta}.treefile"
    mpboot_log = "${fasta}.log"
    tree_w_branchlength = "${pdid}.tree_with_branch_length.tree"
    tree_w_branchlenth_plot = "${pdid}.tree_with_branch_length.pdf"
    phylogeny_vcf = "${pdid}.muts_assigned_to_tree.vcf"
    """
    genotype_bin_to_fasta.py --genotype_bin_path ${genotype_bin_path} --outfile ${fasta}
    mpboot -s $fastafile -bb 1000 -cost e -st DNA
    Rscript --vanilla $script --nr_path=$nr_path --nv_path=$nv_path --genotype_bin_path=$genotype_bin_path --tree_path=$tree_path --out_prefix=$pdid
    """

}