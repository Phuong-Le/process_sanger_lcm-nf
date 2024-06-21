#!/usr/bin/env nextflow 

// Enable DSL 2
nextflow.enable.dsl=2

// import processes for main workflow
include { validate } from "$projectDir/modules/validate.nf"

// include different workflow options
include { WITH_MATCH_NORMAL_SNP } from "$projectDir/workflows/with_match_normal_snps.nf"



// validate parameters
validate(params)
// download container images
include { singularityPreflight } from "$projectDir/modules/singularity"
// If Singularity is used as the container engine and not showing help message, do preflight check to prevent parallel pull issues
// Related issue: https://github.com/nextflow-io/nextflow/issues/1210
if (workflow.containerEngine == 'singularity') {
    singularityPreflight(workflow.container, params.singularity_cachedir)
}


workflow {
    if (params.with_match_normal == true & params.mut_type == 'snp') {
        WITH_MATCH_NORMAL_SNP()
    }

    
}