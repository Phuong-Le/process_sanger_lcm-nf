process cgpVaf {
    publishDir "${params.outdir}/${pdid}", overwrite: false

    // calculate Vaf 
    input:
    tuple val(pdid), val(sample_id_ls), val(match_normal_id_ls), path(vcf_ls), path(tbi_ls), path(bam_ls), path(bai_ls), path(bam_match_ls, stageAs: 'bam_match*.bam'), path(bai_match_ls, stageAs: 'bam_match*.bam.bai')
    val mut_type
    // path reference_genome // use abolute path in cmd below as cgpVaf can only take absolute values // keeping this in case cgpVaf module changes such that absolute path is no longer required
    // path high_depth_region 

    output:
    tuple val(pdid), val(sample_id_ls), val(match_normal_id), path("${pdid}_${mut_type}_vaf.tsv", arity: '1')

    script:
    // assuming the matching normal files are identical, the first file in the list will be used 
    def bam_match = bam_match_ls[0]
    // similarly for match normal id 
    match_normal_id = match_normal_id_ls[0]

    // sorting VCF files in the same order as the sample ID, assumine sample_ID is used as the simple name
    vcf_ls.sort { a, b ->
        indexA = sample_id_ls.indexOf(a.getSimpleName().toString())
        indexB = sample_id_ls.indexOf(b.getSimpleName().toString())
        indexA <=> indexB
    } 
    vcf = vcf_ls.join(" ")

    // sorting BAM files in the same order as the sample ID, assumine sample_ID is used as the simple name
    bam_ls.sort { a, b ->
        def indexA = sample_id_ls.indexOf(a.getSimpleName().toString())
        def indexB = sample_id_ls.indexOf(b.getSimpleName().toString())
        indexA <=> indexB
    } 
    bam = bam_ls.join(" ")

    // creating cgpVaf command
    samples = sample_id_ls.join(" ")
    cmd="cgpVaf.pl -d . -o . -a  ${mut_type} -g  ${params.reference_genome} -hdr  ${params.high_depth_region} --vcf  ${vcf} --normal_bam ${bam_match} --tumour_bam ${bam} --normal_name ${match_normal_id} --tumour_name ${samples} -ct 1" // params.reference_genome and params.high_depth_region because cgpVaf can only take absolute path
    

    """
    $cmd
    grep -vE '^(##)' "${match_normal_id}_${sample_id_ls[0]}_${mut_type}_vaf.tsv" | sed 's/#//' > "${pdid}_${mut_type}_vaf.tsv"
    """
    
}