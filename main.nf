#!/usr/bin/env nextflow 

// Enable DSL 2
nextflow.enable.dsl=2


// import processes for main workflow
include { validate } from "$projectDir/modules/validate.nf"
include { hairpinFilter } from "$projectDir/modules/hairpin.nf"
include { cgpVaf } from "$projectDir/modules/cgpVaf.nf"
include { betaBinomFilterIndex } from "$projectDir/modules/betaBinomFilterIndex.nf"
include { betaBinomFilter } from "$projectDir/modules/betaBinomFilter.nf"

// validate parameters
validate(params)

workflow {
    String sample_paths = new File(params.sample_paths).getText('UTF-8')

    // Hairpin filtering for SNPs 
    if (params.mut_type=='snp') {
        hairpin_input_ch = Channel.of(sample_paths)
            .splitCsv( header: true, sep : '\t' )
            .map { row -> tuple( row.sample_id, row.match_normal_id, row.pdid, row.vcf, row.vcf_tbi, row.bam, row.bai, row.bas, row.met, row.bam_match, row.bai_match ) }
        hairpin_filtered_ch = hairpinFilter(hairpin_input_ch, params.vcfilter_config)

        // group by PDID for cgpVaf
        cgpvaf_input_ch = hairpin_filtered_ch.groupTuple( by: 0 )
    }

    // Quality filtering based on Pindel for Indels 
    
    // cgpVaf 
    cgpVaf_out_ch = cgpVaf(cgpvaf_input_ch, params.mut_type)
    // cgpVaf_out_ch = cgpVaf(cgpvaf_input_ch, params.mut_type, params.reference_genome, params.high_depth_region) // keeping this in case cgpVaf module changes such that absolute path is no longer required

    (beta_binom_index_ch, germline, somatic, rho) = betaBinomFilterIndex(cgpVaf_out_ch)

    // filtering from the hairpin filtered file 
    hairpin_filtered_relevant_ch = hairpin_filtered_ch
        .map( sample -> tuple(sample[0], sample[1], sample[2], sample[3], sample[4]) )

    beta_binom_filter_input_ch = beta_binom_index_ch.cross(hairpin_filtered_relevant_ch)
        .map( sample -> tuple(sample[0][0], sample[1][1], sample[1][2], sample[1][3], sample[1][4], sample[0][1]) )
        
    out = betaBinomFilter(beta_binom_filter_input_ch)
    out.view()

}