# Filter_sanger_lcm-nf

A Nextflow pipeline to process LCM data post-CaVEMan or Pindel variant calling 

## Dependencies

- [Nextflow](https://www.nextflow.io/)
- [hairpin](https://github.com/cancerit/hairpin-wrapper) wrapper by CASM IT at the Sanger institute
- [vcfilter]() - a wrapper by CASM IT (no public repo at the moment but I might switch to bcftools)
- [tabix](https://www.htslib.org/doc/tabix.html)
- [cgpVAFcommand](https://github.com/cancerit/vafCorrect/tree/dev)
- [R](https://www.r-project.org/)

## Installation

Clone this repository

```
git clone git@github.com:Phuong-Le/process_sanger_lcm-nf.git
```

## Usage 

### Arguments

""" 

`sample_paths`:       a file that contains paths to relevant files for the processes in the pipeline, \
    with headers similar to [demo_files/lcm_processing_input_indel.tsv](https://github.com/Phuong-Le/process_sanger_lcm-nf/blob/main/demo_files/lcm_processing_input_indel.tsv) for Indels \
    or [demo_files/lcm_processing_input_snv.tsv](https://github.com/Phuong-Le/process_sanger_lcm-nf/blob/main/demo_files/lcm_processing_input_snv.tsv) for SNVs. \
    NOTE that the `data_dir` column is not neccessary. \
        the [Simple Name](https://www.geeksforgeeks.org/class-getsimplename-method-in-java-with-examples/) of the files should be the same as the sample ID \
        `vcf` filenames need to end with "vcf.gz" \
        Due to cgpVaf requirements, the bam files should be named "${sample_id}.bam" \
        the bai files should be named "${sample_id}.bam.bai" \
        - similarly for the matching normal sample "${match_normal_id}.bam" and "${match_normal_id}.bam.bai"

`vcfilter_config`:    pre-cgpVaf filtering, default to [data/indel_default.filter](https://github.com/Phuong-Le/process_sanger_lcm-nf/blob/main/data/indel_default.filter) for Indels and [data/snv_default.filter](https://github.com/Phuong-Le/process_sanger_lcm-nf/blob/main/data/snv_default.filter) for SNVs

`mut_type`:           either "snv" (default) or "indel"

`reference_genome`:   absolute path to reference genome fasta file

`high_depth_region`:  absolute path to high depth region bed file

`reference_genome_cachedir`: cachedir for reference genome - if this directory is empty the pipeline will split sequences from `reference_genome`, make sure this directory has reference file for all values 'CHROM' column in the produced VCF file
    
`mutmat_kmer`: the kmer size to construct the mutation matrix, default to 3

`outdir`:             path to out directory

"""

### Examples

```
nf_script=path/to/process_sanger_lcm-nf/main.nf
config_file=path/to/config/file
sample_paths=path/to/lcm_processing_input_snv.tsv
mut_type=snv
reference_genome=absolute/path/to/genome.fa
high_depth_region=absolute/path/to/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz
reference_genome_cachedir=path/to/genome/cachedir
mutmat_kmer=3
outdir=path/to/out/directory


nextflow run ${nf_script} -c ${config_file} \
    --sample_paths ${sample_paths} \
    --mut_type ${mut_type} \
    --reference_genome ${reference_genome} \
    --high_depth_region ${high_depth_region} \
    --outdir ${outdir}"
```

## Usage for Sanger's Farm users

### Examples

```
module load nextflow/23.10.1-5891

module load hairpin/1.0.6
module load vcfilter/1.0.4
module load tabix/1.13

module load cgpVAFcommand/2.5.0

module load R/4.1.0 

working_dir=path/to/working/directory

nf_script=path/to/process_sanger_lcm-nf/main.nf
config_file=path/to/process_sanger_lcm-nf/sanger_lsf.config

sample_paths=path/to/lcm_processing_input_snv.tsv
mut_type=snv

reference_genome=/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa
high_depth_region=/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz
reference_genome_cachedir=path/to/genome/cachedir
mutmat_kmer=3
outdir=${working_dir}

bsub -cwd ${working_dir} -o %J.${mut_type}.out -e %J.${mut_type}.err -R "select[mem>1000] rusage[mem=1000]" -M1000 \
    "nextflow run ${nf_script} -c ${config_file} --sample_paths ${sample_paths} --mut_type ${mut_type} --reference_genome ${reference_genome} --high_depth_region ${high_depth_region} --outdir ${outdir}"

```

## Detailed description of the pipeline 

Step 1: Filtering based on CaVeMan or Pindels. (default filtering criteria can be found in either [data/indel_default.filter](https://github.com/Phuong-Le/process_sanger_lcm-nf/blob/main/data/indel_default.filter) or [data/snv_default.filter](https://github.com/Phuong-Le/process_sanger_lcm-nf/blob/main/data/snv_default.filter))

Step 2: for each donor with multiple samples, use cgpVaf to calculate VAF across their samples.

Step 3: Use beta-binomial test based on VAF to filter out germline and LCM artefact mutations.

Step 4: Generate mutation matrix (only implemented for snps at the moment, will add this feature for indels)

## Authors 

[Phuong Le](https://github.com/Phuong-Le) (email: al35@sanger.ac.uk)

## Notes

This pipeline (particularly the sanger_lsf.config) has only be tested on non-malignant colon data (usually 2000-9000 mutations per sample). If you have a bigger dataset, it might be neccessary to experiment with different memory settings. 

If neccessary, future work might \
- improve cgpVaf by parallelising computation on chromosome


## Acknowledgements

Thanks to Chloe Pacyna and [Yichen Wang] for guiding me through the necessary steps to process LCM data

Thanks to Shriram Bhosle from CASM IT for showing me how to use cgpVaf (ie vafCorrect)

The BetaBinomial filtering of germline and artefacts was adapted from Tim Coorens's code at https://github.com/TimCoorens/Unmatched_NormSeq/tree/master

The validatation of parameter was adapted from code by Harry Hung