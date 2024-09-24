#!/usr/bin/env nextflow 

// enable DSL 2
nextflow.enable.dsl=2

// all of the default parameters are set in `nextflow.config`

// include functions
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// include different workflow options
include { CONPAIR_FILTER_WITH_MATCH_NORMAL } from "$projectDir/workflows/conpair_filter_with_match_normal.nf"
include { FILTER_WITH_MATCH_NORMAL_SNV } from "$projectDir/workflows/filter_with_match_normal_snvs.nf"
include { FILTER_WITH_MATCH_NORMAL_INDEL } from "$projectDir/workflows/filter_with_match_normal_indels.nf"
include { PHYLOGENETICS } from "$projectDir/workflows/phylogenetics.nf"
include { PHYLOGENETICS_PROVIDED_TREE_TOPOLOGY } from "$projectDir/workflows/phylogenetics_provided_topology.nf"

// download container images
include { singularityPreflight } from "$projectDir/modules/singularity"
// If Singularity is used as the container engine and not showing help message, do preflight check to prevent parallel pull issues
// Related issue: https://github.com/nextflow-io/nextflow/issues/1210
if (workflow.containerEngine == 'singularity') {
    singularityPreflight(workflow.container, params.singularity_cachedir)
}

// print help message, with typical command line usage
if (params.help) {
  def String command = """nextflow run process_sanger_lcm-nf \\
    --samplesheet /path/to/samplesheet.csv \\
    --reference_genome /path/to/genome.fa \\
    --high_depth_bed /path/to/HiDepth.bed.gz \\
    --outdir out/""".stripIndent()
  log.info paramsHelp(command)
  exit 0
}

// validate parameters
validateParameters()

workflow {
    
    // print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    // generate a channel of samples
    //ch_input = Channel.fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
    Channel.fromPath(params.samplesheet, checkIfExists: true)
    | splitCsv(header: true)
    | map { row ->
            meta = row.subMap("sample_id", "match_normal_id", "pdid")
            [meta,
              file(row.bam),
              file(row.bai),
              file(row.bas),
              file(row.met),
              file(row.vcf_snv),
              file(row.vcf_tbi_snv),
              file(row.vcf_indel),
              file(row.vcf_tbi_indel),
              file(row.bam_match),
              file(row.bai_match)]
    }
    | set { ch_input }
    ch_input.view()
    exit 0

    println("test")
    exit 0

    if (params.with_match_normal == true) {

        // conpair 
        if (params.conpair == true) {
            CONPAIR_FILTER_WITH_MATCH_NORMAL(params.samplesheet) 
        }

        // filtering snvs 
        if (params.filter_snv == true) {
            vcfilter_config = (params.vcfilter_config=="") ? "${projectDir}/data/snv_default.filter" : params.vcfilter_config
            if (params.conpair == true) {
                samplesheet_content_ch = CONPAIR_FILTER_WITH_MATCH_NORMAL.out
                    .splitCsv( header: true, sep : '\t' )
                    .map { row -> tuple( row.sample_id, row.match_normal_id, row.pdid, row.vcf_snv, row.vcf_tbi_snv, row.bam, row.bai, row.bas, row.met, row.bam_match, row.bai_match ) }
                FILTER_WITH_MATCH_NORMAL_SNV(samplesheet_content_ch, vcfilter_config)
            }
            else {
                samplesheet = new File(params.samplesheet).getText('UTF-8')
                samplesheet_content_ch = Channel.of(samplesheet)
                    .splitCsv( header: true, sep : '\t' )
                    .map { row -> tuple( row.sample_id, row.match_normal_id, row.pdid, row.vcf_snv, row.vcf_tbi_snv, row.bam, row.bai, row.bas, row.met, row.bam_match, row.bai_match ) }
                FILTER_WITH_MATCH_NORMAL_SNV(samplesheet_content_ch, vcfilter_config)
            }
        }

        // filtering indels 
        if (params.filter_indel == true) {
            vcfilter_config = (params.vcfilter_config=="") ? "${projectDir}/data/indel_default.filter" : params.vcfilter_config
            if (params.conpair == true) {
                samplesheet_content_ch = CONPAIR_FILTER_WITH_MATCH_NORMAL.out
                    .splitCsv( header: true, sep : '\t' )
                    .map { row -> tuple( row.sample_id, row.match_normal_id, row.pdid, row.vcf_indel, row.bam, row.bai, row.bam_match, row.bai_match ) }
                FILTER_WITH_MATCH_NORMAL_INDEL(samplesheet_content_ch, vcfilter_config)
            }
            else {
                samplesheet = new File(params.samplesheet).getText('UTF-8')
                samplesheet_content_ch = Channel.of(samplesheet)
                    .splitCsv( header: true, sep : '\t' )
                    .map { row -> tuple( row.sample_id, row.match_normal_id, row.pdid, row.vcf_indel, row.bam, row.bai, row.bam_match, row.bai_match ) }
                FILTER_WITH_MATCH_NORMAL_INDEL(samplesheet_content_ch, vcfilter_config)
            }
        }
        
    }

    // phylogenetics is independent of whether there's a match normal or not
    // indel phylogenetics will use output from snv phylogenetics if both workflows are run
    if (params.phylogenetics == true) {
        if (params.filter_snv == true) {
            // only run this if there are more than 2 sample per donor (genotype_bin only has one column)
            phylogenetics_input_ch = FILTER_WITH_MATCH_NORMAL_SNV
                .out
                .filter { it[3].readLines().first().split(' ').size() > 2 }
            PHYLOGENETICS(phylogenetics_input_ch, 'phylogenetics_snv_out') // phylogenetics without tree topology
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
            samplesheet = new File(params.samplesheet).getText('UTF-8')
            topology = Channel.of(samplesheet)
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
            assert params.with_topology != null 
            if (params.with_topology == true) {
                // process input samplesheet
                outdir_basename = (params.phylogenetics_outdir_basename == "") ? 'phylogenetics_indel_out' : params.phylogenetics_outdir_basename
                samplesheet = new File(params.samplesheet).getText('UTF-8')
                sample_path_content = Channel.of(samplesheet)
                    .splitCsv( header: true, sep : '\t' )
                    .map{ row -> tuple( row.pdid, row.topology, row.nr_path, row.nv_path, row.genotype_bin_path ) }
                PHYLOGENETICS_PROVIDED_TREE_TOPOLOGY(sample_path_content, outdir_basename)
            }
            else {
                // process input sample paths 
                outdir_basename = (params.phylogenetics_outdir_basename == "") ? 'phylogenetics_snv_out' : params.phylogenetics_outdir_basename
                samplesheet = new File(params.samplesheet).getText('UTF-8')
                sample_path_content = Channel.of(samplesheet)
                    .splitCsv( header: true, sep : '\t' )
                    .map { row -> tuple( row.pdid, row.nr_path, row.nv_path, row.genotype_bin_path ) }
                PHYLOGENETICS(sample_path_content, outdir_basename)
            }
        }
    }
       
}