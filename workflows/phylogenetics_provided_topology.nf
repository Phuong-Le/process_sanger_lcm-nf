include { mutToTree } from "$projectDir/modules/mutToTree.nf"
include { matrixGeneratorOnBranches } from "$projectDir/modules/matrixGeneratorOnBranches.nf"
include { concatMatrices } from "$projectDir/modules/concatMatrices.nf"
include { sigprofilerPlotSnpByBranches } from "$projectDir/modules/sigprofilerPlotSnpByBranches.nf"

workflow PHYLOGENETICS_PROVIDED_TREE_TOPOLOGY { // phylogenetics workflow for INDELs
    take:
    mutToTree_input_ch
    outdir_basename
    
    main:
    // assign mutation to provided topology
    (branched_vcf_with_header, other_files) = mutToTree(mutToTree_input_ch, outdir_basename)

    // generate mutation matrix for the branches by SigProfilerMatrixGenerator
    matrixGeneratorOnBranches(branched_vcf_with_header, outdir_basename)
    concatMatrices(matrixGeneratorOnBranches.out.toList(), outdir_basename)
    // plotting
    // sigprofilerPlotSnpByBranches(concatMatrices.out, outdir_basename)

}