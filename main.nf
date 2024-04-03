#!/usr/bin/env nextflow 

// Enable DSL 2
nextflow.enable.dsl=2

// import processes for main workflow
include { validate } from "$projectDir/modules/validate.nf"
include { conpairPileup as conpairPileupSample } from "$projectDir/modules/conpairPileup.nf"
include { conpairPileup as conpairPileupMatch } from "$projectDir/modules/conpairPileup.nf"
include { verifyConcordance } from "$projectDir/modules/verifyConcordance.nf"
include { conpairContamination } from "$projectDir/modules/conpairContamination.nf"
include { hairpinFilter } from "$projectDir/modules/hairpin.nf"
include { pindelFilter } from "$projectDir/modules/pindelFilter.nf"
include { cgpVaf } from "$projectDir/modules/cgpVaf.nf"
include { betaBinomFilterIndex } from "$projectDir/modules/betaBinomFilterIndex.nf"
include { betaBinomFilter } from "$projectDir/modules/betaBinomFilter.nf"
include { splitRef } from "$projectDir/modules/splitRef.nf"
include { getMutMat } from "$projectDir/modules/getMutMat.nf"

// validate parameters
validate(params)

// if reference genome cachedir doesn't exist then create one 
if ( !file(params.reference_genome_cachedir).exists() ) { 
    params.reference_genome_cachedir.mkdir()
    }  

workflow {
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
    // concordance_input_ch.view()
    concordance_output_ch = verifyConcordance(concordance_input_ch)
        .collectFile( name: 'concordance.txt', newLine: true )
    concordance_output_ch.view()

    // Contamination 
    contamination_input_ch = pileup_sample.cross(pileup_match) 
        .map { sample -> tuple(sample[0][0], sample[1][2], sample[0][1], sample[0][2] ) }
    contamination_input_ch.view()
    contamination_output_ch = conpairContamination(contamination_input_ch)
        .collectFile( name: 'contamination.txt', newLine: true )
    contamination_output_ch.view()

    // Hairpin filtering for SNPs 
    // if (params.mut_type=='snp') {
    //     vcfilter_config = (params.vcfilter_config=="") ? "${projectDir}/data/snp_default.filter" : params.vcfilter_config
    //     hairpin_input_ch = Channel.of(sample_paths)
    //         .splitCsv( header: true, sep : '\t' )
    //         .map { row -> tuple( row.sample_id, row.match_normal_id, row.pdid, row.vcf, row.vcf_tbi, row.bam, row.bai, row.bas, row.met, row.bam_match, row.bai_match ) }
    //     vcfiltered_ch = hairpinFilter(hairpin_input_ch, vcfilter_config)
    // }
    
    // // // Quality filtering based on Pindel for Indels 
    // // if (params.mut_type=='indel') {
    // //     vcfilter_config = (params.vcfilter_config=="") ? "${projectDir}/data/indel_default.filter" : params.vcfilter_config
    // //     indel_input_ch = Channel.of(sample_paths)
    // //         .splitCsv( header: true, sep : '\t' )
    // //         .map { row -> tuple( row.sample_id, row.match_normal_id, row.pdid, row.vcf, row.bam, row.bai, row.bam_match, row.bai_match ) }
    // //     vcfiltered_ch = pindelFilter(indel_input_ch, vcfilter_config)
    // // }


    // // cgpVaf 
    // cgpvaf_input_ch = vcfiltered_ch.groupTuple( by: 0 ) // group by PDID for cgpVaf
    // cgpVaf_out_ch = cgpVaf(cgpvaf_input_ch, params.mut_type)
    // // cgpVaf_out_ch = cgpVaf(cgpvaf_input_ch, params.mut_type, params.reference_genome, params.high_depth_region) // keeping this in case cgpVaf module changes such that absolute path is no longer required

    // // BetaBinomial filtering for germline and LCM artefacts based on cgpVaf (methods by Tim Coorens)
    // (beta_binom_index_ch, germline, somatic, rho) = betaBinomFilterIndex(cgpVaf_out_ch) // get the indices for the filtering 
    // // use hairpin or pindel vcfiltered output to recover the donor-based channels from cgpVaf
    // vcfiltered_relevant_ch = vcfiltered_ch
    //     .map( sample -> tuple(sample[0], sample[1], sample[2], sample[3], sample[4]) )
    // beta_binom_filter_input_ch = beta_binom_index_ch.cross(vcfiltered_relevant_ch)
    //     .map( sample -> tuple(sample[0][0], sample[1][1], sample[1][2], sample[1][3], sample[1][4], sample[0][1]) )
    // (bbinom_filtered_vcf_ch, filtered_vcf_noheader_ch) = betaBinomFilter(beta_binom_filter_input_ch)
    
    // // split reference genome if not cached (ie if cachedir is empty)
    // if ( file(params.reference_genome_cachedir).listFiles().toList().isEmpty() ) {
    //     ref_chrom_ch = Channel.fromPath(params.reference_genome)
    //         .splitFasta( file: true )
    //     split_fastas = splitRef(ref_chrom_ch).toList()
    // } else { 
    //     split_fastas = file(params.reference_genome_cachedir).listFiles().toList() 
    // }
    // // only the main chromosomes are allowed
    // chromosomes = (1..22).collect { "chr${it}".toString() } + ['chrX', 'chrY']
    // split_fastas_chrom = split_fastas.findAll { 
    //     chromosomes.contains("${it}".tokenize("/").last().tokenize(".").first() ) 
    //     }
    
    // // generate the mutation matrix - not yet available for Indels
    // if (params.mut_type == 'snp') {
    //     mutation_matrix_ch = getMutMat(bbinom_filtered_vcf_ch, split_fastas_chrom, params.mutmat_kmer)
    //     mutation_matrix_ch.collectFile(name: "${params.outdir}/mutation_matrix.txt", keepHeader: true, skip: 1)
    // }
}