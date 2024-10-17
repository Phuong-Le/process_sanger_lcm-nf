#!/usr/bin/env nextflow 

// Enable DSL 2
nextflow.enable.dsl=2

// import processes for main workflow
include { validate } from "$projectDir/modules/validate.nf"

// include different workflow options
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { CONPAIR_FILTER_WITH_MATCH_NORMAL } from "$projectDir/workflows/conpair_filter_with_match_normal.nf"
include { FILTER_WITH_MATCH_NORMAL_SNP } from "$projectDir/workflows/filter_with_match_normal_snps.nf"
include { FILTER_WITH_MATCH_NORMAL_INDEL } from "$projectDir/workflows/filter_with_match_normal_indels.nf"
include { PHYLOGENETICS } from "$projectDir/workflows/phylogenetics.nf"
include { PHYLOGENETICS_PROVIDED_TREE_TOPOLOGY } from "$projectDir/workflows/phylogenetics_provided_topology.nf"

// validate parameters
// validate(params)
// download container images
include { singularityPreflight } from "$projectDir/modules/singularity"
// If Singularity is used as the container engine and not showing help message, do preflight check to prevent parallel pull issues
// Related issue: https://github.com/nextflow-io/nextflow/issues/1210
if (workflow.containerEngine == 'singularity') {
    singularityPreflight(workflow.container, params.singularity_cachedir)
}

// validate parameters
validateParameters()

workflow {

    // String sample_paths = new File(params.sample_paths).getText('UTF-8')
    
    if (params.with_match_normal == true) {

        // conpair 
        if (params.conpair == true) {
            CONPAIR_FILTER_WITH_MATCH_NORMAL(params.sample_paths) 
        }

        // filtering snps 
        if (params.filter_snp == true) {
            if (params.conpair == true) {
                sample_paths_content_ch = CONPAIR_FILTER_WITH_MATCH_NORMAL.out
                    .splitCsv( header: true, sep : '\t' )
                    .map { row -> tuple( row.sample_id, row.match_normal_id, row.pdid, row.vcf_snp, row.vcf_tbi_snp, row.bam, row.bai, row.bas, row.met, row.bam_match, row.bai_match ) }
                FILTER_WITH_MATCH_NORMAL_SNP(sample_paths_content_ch, params.vcfilter_snv, params.bbinom_rho_snv)
            }
            else {
                sample_paths = new File(params.sample_paths).getText('UTF-8')
                sample_paths_content_ch = Channel.of(sample_paths)
                    .splitCsv( header: true, sep : '\t' )
                    .map { row -> tuple( row.sample_id, row.match_normal_id, row.pdid, row.vcf_snp, row.vcf_tbi_snp, row.bam, row.bai, row.bas, row.met, row.bam_match, row.bai_match ) }
                FILTER_WITH_MATCH_NORMAL_SNP(sample_paths_content_ch, params.vcfilter_snv, params.bbinom_rho_snv)
            }
        }

        // // filtering indels 
        if (params.filter_indel == true) {
            if (params.conpair == true) {
                sample_paths_content_ch = CONPAIR_FILTER_WITH_MATCH_NORMAL.out
                    .splitCsv( header: true, sep : '\t' )
                    .map { row -> tuple( row.sample_id, row.match_normal_id, row.pdid, row.vcf_indel, row.bam, row.bai, row.bam_match, row.bai_match ) }
                FILTER_WITH_MATCH_NORMAL_INDEL(sample_paths_content_ch, params.vcfilter_indel, params.bbinom_rho_indel)
            }
            else {
                sample_paths = new File(params.sample_paths).getText('UTF-8')
                sample_paths_content_ch = Channel.of(sample_paths)
                    .splitCsv( header: true, sep : '\t' )
                    .map { row -> tuple( row.sample_id, row.match_normal_id, row.pdid, row.vcf_indel, row.bam, row.bai, row.bam_match, row.bai_match ) }
                FILTER_WITH_MATCH_NORMAL_INDEL(sample_paths_content_ch, params.vcfilter_indel, params.bbinom_rho_indel)
            }
        }

        // Phylogenetics
        if (params.phylogenetics == true) {
            if (params.filter_snp == true) {
                // only run this if there are more than 2 sample per donor (genotype_bin only has one column)
                phylogenetics_input_ch = FILTER_WITH_MATCH_NORMAL_SNP
                    .out
                    .filter { it[3].readLines().first().split(' ').size() > 2 }
                PHYLOGENETICS(phylogenetics_input_ch, 'phylogenetics_snp_out') // phylogenetics without tree topology
                if (params.filter_indel == true) {
                    phylogenetics_input_ch = FILTER_WITH_MATCH_NORMAL_INDEL
                        .out
                        .filter { it[3].readLines().first().split(' ').size() > 2 }
                    mutToTree_input_ch = PHYLOGENETICS.out.cross(phylogenetics_input_ch)
                        .map( pdid -> tuple(pdid[0][0], pdid[0][1], pdid[1][1], pdid[1][2], pdid[1][3]) )
                    PHYLOGENETICS_PROVIDED_TREE_TOPOLOGY(mutToTree_input_ch, 'phylogenetics_indel_out')
                }
            }
            else if (params.filter_indel == true) {
                // get topology
                sample_paths = new File(params.sample_paths).getText('UTF-8')
                topology = Channel.of(sample_paths)
                    .splitCsv( header: true, sep : '\t' )
                    .map { row -> tuple( row.pdid, row.topology ) }
                    .unique()
            
                phylogenetics_input_ch = FILTER_WITH_MATCH_NORMAL_INDEL
                        .out
                        .filter { it[3].readLines().first().split(' ').size() > 2 }

                mutToTree_input_ch = topology.cross(phylogenetics_input_ch)
                        .map( pdid -> tuple(pdid[0][0], pdid[0][1], pdid[1][1], pdid[1][2], pdid[1][3]) )

                PHYLOGENETICS_PROVIDED_TREE_TOPOLOGY(mutToTree_input_ch, 'phylogenetics_indel_out') // phylogenetics with tree topology
            }
            else { // phylogenetics pipeline only 
                if (params.snv_then_indel == true) {
                    sample_paths = new File(params.sample_paths).getText('UTF-8')
                    phylogenetics_snv_input_ch = Channel.of(sample_paths)
                        .splitCsv( header: true, sep : '\t' )
                        .map { row -> tuple( row.pdid, row.nr_path_snv, row.nv_path_snv, row.genotype_bin_path_snv ) }
                    PHYLOGENETICS(phylogenetics_snv_input_ch, 'phylogenetics_snp_out') // phylogenetics for SNV
                    phylogenetics_indel_input_ch = Channel.of(sample_paths)
                        .splitCsv( header: true, sep : '\t' )
                        .map { row -> tuple( row.pdid, row.nr_path_indel, row.nv_path_indel, row.genotype_bin_path_indel ) }
                    mutToTree_input_ch = PHYLOGENETICS.out.cross(phylogenetics_indel_input_ch)
                        .map( pdid -> tuple(pdid[0][0], pdid[0][1], pdid[1][1], pdid[1][2], pdid[1][3]) ) // added topology
                    PHYLOGENETICS_PROVIDED_TREE_TOPOLOGY(mutToTree_input_ch, 'phylogenetics_indel_out')
                }
                else {
                    assert params.with_topology != null
                    if (params.with_topology == true) {
                    // process input sample_paths
                    outdir_basename = (params.phylogenetics_outdir_basename == "") ? 'phylogenetics_indel_out' : params.phylogenetics_outdir_basename
                    sample_paths = new File(params.sample_paths).getText('UTF-8')
                    sample_path_content = Channel.of(sample_paths)
                        .splitCsv( header: true, sep : '\t' )
                        .map{ row -> tuple( row.pdid, row.topology, row.nr_path, row.nv_path, row.genotype_bin_path ) }
                    PHYLOGENETICS_PROVIDED_TREE_TOPOLOGY(sample_path_content, outdir_basename)
                    }
                    else {
                        // process input sample paths 
                        outdir_basename = (params.phylogenetics_outdir_basename == "") ? 'phylogenetics_snp_out' : params.phylogenetics_outdir_basename
                        sample_paths = new File(params.sample_paths).getText('UTF-8')
                        sample_path_content = Channel.of(sample_paths)
                            .splitCsv( header: true, sep : '\t' )
                            .map { row -> tuple( row.pdid, row.nr_path, row.nv_path, row.genotype_bin_path ) }
                        PHYLOGENETICS(sample_path_content, outdir_basename)
                    }
                }
            }
        }
    
    }
       
}