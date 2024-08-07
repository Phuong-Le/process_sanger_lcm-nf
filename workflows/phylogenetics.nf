include { getPhylogeny } from "$projectDir/modules/getPhylogeny.nf"
include { matrixGeneratorOnBranches } from "$projectDir/modules/matrixGeneratorOnBranches.nf"
include { concatMatrices } from "$projectDir/modules/concatMatrices.nf"
include { sigprofilerPlotSnpByBranches } from "$projectDir/modules/sigprofilerPlotSnpByBranches.nf"


workflow PHYLOGENETICS { // phylogenetics workflow for SNVs
    take: 
    phylogenetics_input_ch
    outdir_basename

    main: 
    (branched_vcf_with_header, topology, other_files, mpboot_log) = getPhylogeny(phylogenetics_input_ch, outdir_basename)

    // generate mutation matrix for the branches by SigProfilerMatrixGenerator
    matrixGeneratorOnBranches(branched_vcf_with_header, outdir_basename)
    concatMatrices(matrixGeneratorOnBranches.out.toList(), outdir_basename)
    // plotting
    sigprofilerPlotSnpByBranches(concatMatrices.out, outdir_basename)

    emit:
    topology

}