include { pindelFilter } from "$projectDir/modules/pindelFilter.nf"
include { cgpVaf } from "$projectDir/modules/cgpVaf.nf"
include { betaBinomFilterIndex } from "$projectDir/modules/betaBinomFilterIndex.nf"
include { betaBinomFilter } from "$projectDir/modules/betaBinomFilter.nf"
include { matrixGeneratorOnSamples } from "$projectDir/modules/matrixGeneratorOnSamples.nf"
include { sigprofilerPlotSnvBySamples } from "$projectDir/modules/sigprofilerPlotSnvBySamples.nf"



workflow FILTER_WITH_MATCH_NORMAL_INDEL {
    take:
    samplesheet_content_ch
    vcfilter_config

    main:
    // setup
    mut_type = 'indel'
    

    // FILTER
    vcfiltered_ch = pindelFilter(samplesheet_content_ch, vcfilter_config, mut_type)


    // cgpVaf 
    cgpVaf_out_ch = cgpVaf(vcfiltered_ch.groupTuple( by: 0 ), mut_type)
    // cgpVaf_out_ch = cgpVaf(cgpvaf_input_ch, params.mut_type, params.reference_genome, params.high_depth_bed) // keeping this in case cgpVaf module changes such that absolute path is no longer required

    // BetaBinomial filtering for germline and LCM artefacts based on cgpVaf (methods by Tim Coorens)
    (beta_binom_index_ch, germline, somatic, rho, phylogenetics_input_ch) = betaBinomFilterIndex(cgpVaf.out, mut_type) // get the indices for the filtering 
    // use hairpin vcfiltered output to recover the donor-based channels from cgpVaf
    vcfiltered_relevant_ch = vcfiltered_ch
        .map( sample -> tuple(sample[0], sample[1], sample[2], sample[3], sample[4]) )
    beta_binom_filter_input_ch = beta_binom_index_ch.cross(vcfiltered_relevant_ch)
        .map( sample -> tuple(sample[0][0], sample[1][1], sample[1][2], sample[1][3], sample[1][4], sample[0][1]) )
    betaBinomFilter(beta_binom_filter_input_ch, mut_type)

    // generate mutation matrix for the samples by SigProfilerMatrixGenerator
    matrixGeneratorOnSamples(betaBinomFilter.out.toList(), mut_type)

    // plot spectra
    sigprofilerPlotSnvBySamples(matrixGeneratorOnSamples.out, mut_type)

    emit:
    phylogenetics_input_ch 
}