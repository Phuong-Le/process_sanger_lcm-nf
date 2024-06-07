include { conpairPileup as conpairPileupSample } from "$projectDir/modules/conpairPileup.nf"
include { conpairPileup as conpairPileupMatch } from "$projectDir/modules/conpairPileup.nf"
include { conpairFilter } from "$projectDir/modules/conpairFilter.nf"
include { verifyConcordance } from "$projectDir/modules/verifyConcordance.nf"
include { conpairContamination } from "$projectDir/modules/conpairContamination.nf"
include { hairpinFilter } from "$projectDir/modules/hairpin.nf"
include { pindelFilter } from "$projectDir/modules/pindelFilter.nf"
include { cgpVaf } from "$projectDir/modules/cgpVaf.nf"
include { betaBinomFilterIndex } from "$projectDir/modules/betaBinomFilterIndex.nf"
include { betaBinomFilter } from "$projectDir/modules/betaBinomFilter.nf"
include { getPhylogeny } from "$projectDir/modules/getPhylogeny.nf"
include { matrixGeneratorOnSamples } from "$projectDir/modules/matrixGeneratorOnSamples.nf"
include { matrixGeneratorOnBranches } from "$projectDir/modules/matrixGeneratorOnBranches.nf"
include { concatMatrices } from "$projectDir/modules/concatMatrices.nf"


workflow WITH_MATCH_NORMAL_SNP {
    main:
    String sample_paths = new File(params.sample_paths).getText('UTF-8')

    // Conpair 
    // sample
    sample_pileup_input_ch = Channel.of(sample_paths)
            .splitCsv( header: true, sep : '\t' )
            .map { row -> tuple( row.match_normal_id, row.sample_id, row.bam, row.bai ) }
    pileup_sample = conpairPileupSample(sample_pileup_input_ch)
    // normal
    match_pileup_input_ch = Channel.of(sample_paths)
            .splitCsv( header: true, sep : '\t' )
            .map { row -> tuple( row.match_normal_id, row.match_normal_id, row.bam_match, row.bai_match ) }
            .unique()
    pileup_match = conpairPileupMatch(match_pileup_input_ch)

    // Concordance between sample and match normal
    concordance_input_ch = pileup_sample.combine(pileup_match)
        .map { sample -> tuple(sample[1], sample[2], sample[3], sample[5]) }
    concordance_output_ch = verifyConcordance(concordance_input_ch)
        .collectFile( name: 'concordance.txt', newLine: true )

    // Contamination 
    contamination_input_ch = pileup_match.cross(pileup_sample)
        .map { sample -> tuple(sample[0][0], sample[0][2], sample[1][1], sample[1][2]) }
    contamination_output_ch = conpairContamination(contamination_input_ch)
        .collectFile( name: 'contamination.txt', newLine: true )

    // Filtering contamination based on concordance and contamination
    (sample_paths_conpaired, conpair_log, concordance_path, contamination_path) = conpairFilter(concordance_output_ch, contamination_output_ch, params.sample_paths)


    // filtered sample paths will replace sample paths from here
    
    // Hairpin filtering for SNPs 
    if (params.mut_type=='snp') {
        vcfilter_config = (params.vcfilter_config=="") ? "${projectDir}/data/snp_default.filter" : params.vcfilter_config
        hairpin_input_ch = sample_paths_conpaired
            .splitCsv( header: true, sep : '\t' )
            .map { row -> tuple( row.sample_id, row.match_normal_id, row.pdid, row.vcf, row.vcf_tbi, row.bam, row.bai, row.bas, row.met, row.bam_match, row.bai_match ) }
        vcfiltered_ch = hairpinFilter(hairpin_input_ch, vcfilter_config)
    }
    
    // Quality filtering based on Pindel for Indels 
    if (params.mut_type=='indel') {
        vcfilter_config = (params.vcfilter_config=="") ? "${projectDir}/data/indel_default.filter" : params.vcfilter_config
        indel_input_ch = sample_paths_conpaired
            .splitCsv( header: true, sep : '\t' )
            .map { row -> tuple( row.sample_id, row.match_normal_id, row.pdid, row.vcf, row.bam, row.bai, row.bam_match, row.bai_match ) }
        vcfiltered_ch = pindelFilter(indel_input_ch, vcfilter_config)
    }


    // cgpVaf 
    cgpvaf_input_ch = vcfiltered_ch.groupTuple( by: 0 ) // group by PDID for cgpVaf
    cgpVaf_out_ch = cgpVaf(cgpvaf_input_ch, params.mut_type)
    // cgpVaf_out_ch = cgpVaf(cgpvaf_input_ch, params.mut_type, params.reference_genome, params.high_depth_region) // keeping this in case cgpVaf module changes such that absolute path is no longer required

    // BetaBinomial filtering for germline and LCM artefacts based on cgpVaf (methods by Tim Coorens)
    (beta_binom_index_ch, germline, somatic, rho, phylogenetics_input_ch) = betaBinomFilterIndex(cgpVaf_out_ch) // get the indices for the filtering 
    // use hairpin or pindel vcfiltered output to recover the donor-based channels from cgpVaf
    vcfiltered_relevant_ch = vcfiltered_ch
        .map( sample -> tuple(sample[0], sample[1], sample[2], sample[3], sample[4]) )
    beta_binom_filter_input_ch = beta_binom_index_ch.cross(vcfiltered_relevant_ch)
        .map( sample -> tuple(sample[0][0], sample[1][1], sample[1][2], sample[1][3], sample[1][4], sample[0][1]) )
    (bbinom_filtered_vcf_ch, filtered_sigprofiler_vcf_ch) = betaBinomFilter(beta_binom_filter_input_ch)


    // Phylogenetics, only run this if there are more than 2 sample per donor (genotype_bin only has one column), AND if mutation type is snp
    // work in progress
    if (params.phylogenetics == true) {
        phylogenetics_input__filtered_ch = phylogenetics_input_ch.filter { it[3].readLines().first().split(' ').size() > 2 }
        phylogenetics_input__filtered_ch.view()
        (branched_vcf, other_files, mpboot_log) = getPhylogeny(phylogenetics_input_ch)
        // generate mutation matrix for the branches by SigProfilerMatrixGenerator
        (matrix_by_branches_ch, vcf_with_header_ch) = matrixGeneratorOnBranches(branched_vcf)
        matrix_by_branches_ch.view()
        concatMatrices(matrix_by_branches_ch.toList())
    }

    // generate mutation matrix for the samples by SigProfilerMatrixGenerator
    matrixGeneratorOnSamples(filtered_sigprofiler_vcf_ch.toList())
}