#!/usr/bin/env nextflow 

// Enable DSL 2
nextflow.enable.dsl=2

// import processes for main workflow
include { validate } from "$projectDir/modules/validate.nf"

// include different workflow options
include { WITH_MATCH_NORMAL_SNP } from "$projectDir/workflows/with_match_normal_snps.nf"



// include { conpairPileup as conpairPileupSample } from "$projectDir/modules/conpairPileup.nf"
// include { conpairPileup as conpairPileupMatch } from "$projectDir/modules/conpairPileup.nf"
// include { conpairFilter } from "$projectDir/modules/conpairFilter.nf"
// include { verifyConcordance } from "$projectDir/modules/verifyConcordance.nf"
// include { conpairContamination } from "$projectDir/modules/conpairContamination.nf"
// include { hairpinFilter } from "$projectDir/modules/hairpin.nf"
// include { pindelFilter } from "$projectDir/modules/pindelFilter.nf"
// include { cgpVaf } from "$projectDir/modules/cgpVaf.nf"
// include { betaBinomFilterIndex } from "$projectDir/modules/betaBinomFilterIndex.nf"
// include { betaBinomFilter } from "$projectDir/modules/betaBinomFilter.nf"
// include { getPhylogeny } from "$projectDir/modules/getPhylogeny.nf"
// include { matrixGeneratorOnSamples } from "$projectDir/modules/matrixGeneratorOnSamples.nf"
// include { matrixGeneratorOnBranches } from "$projectDir/modules/matrixGeneratorOnBranches.nf"
// include { concatMatrices } from "$projectDir/modules/concatMatrices.nf"
// include { splitRef } from "$projectDir/modules/splitRef.nf"
// include { getMutMat } from "$projectDir/modules/getMutMat.nf"

// validate parameters
validate(params)
// download container images
include { singularityPreflight } from "$projectDir/modules/singularity"
// If Singularity is used as the container engine and not showing help message, do preflight check to prevent parallel pull issues
// Related issue: https://github.com/nextflow-io/nextflow/issues/1210
// if (workflow.containerEngine == 'singularity') {
//     singularityPreflight(workflow.container, params.singularity_cachedir)
// }

// if reference genome cachedir doesn't exist then create one 
// if ( !file(params.reference_genome_cachedir).exists() ) { 
//     params.reference_genome_cachedir.mkdir()
//     }  

workflow {
    if (params.with_match_normal == true & params.mut_type == 'snp') {
        WITH_MATCH_NORMAL_SNP()
    }

    

//     // // split reference genome if not cached (ie if cachedir is empty)
//     // if ( file(params.reference_genome_cachedir).listFiles().toList().isEmpty() ) {
//     //     ref_chrom_ch = Channel.fromPath(params.reference_genome)
//     //         .splitFasta( file: true )
//     //     split_fastas = splitRef(ref_chrom_ch).toList()
//     // } else { 
//     //     split_fastas = file(params.reference_genome_cachedir).listFiles().toList() 
//     // }
//     // // only the main chromosomes are allowed
//     // chromosomes = (1..22).collect { "chr${it}".toString() } + ['chrX', 'chrY']
//     // split_fastas_chrom = split_fastas.findAll { 
//     //     chromosomes.contains("${it}".tokenize("/").last().tokenize(".").first() ) 
//     //     }
    
//     // // generate the mutation matrix - not yet available for Indels
//     // if (params.mut_type == 'snp') {
//     //     mutation_matrix_ch = getMutMat(bbinom_filtered_vcf_ch, split_fastas_chrom, params.mutmat_kmer)
//     //     mutation_matrix_ch.collectFile(name: "${params.outdir}/mutation_matrix.txt", keepHeader: true, skip: 1)
//     // }
}