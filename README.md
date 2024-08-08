# process_sanger_lcm-nf

## Introduction

`process_sanger_lcm-nf` is a [Nextflow](https://www.nextflow.io/) pipeline to process low-input DNA sequencing data, either from laser capture microdissection or from single-cell derived colonies. It takes as input variant calls in the form of VCF files from `CaVEMAN` (SNVs) or `Pindel` (indels).

1. Filtering based on CaVeMan or Pindels. Default filtering criteria can be found in either [data/indel_default.filter](https://github.com/Phuong-Le/process_sanger_lcm-nf/blob/main/data/indel_default.filter) or [data/snv_default.filter](https://github.com/Phuong-Le/process_sanger_lcm-nf/blob/main/data/snv_default.filter).
2. For each donor with multiple samples, use `cgpVAF` to calculate VAF across their samples.
3. Use beta-binomial test based on VAF to filter out germline and LCM artefact mutations.
4. Generate mutation matrix (only implemented for SNVs at the moment, will add this feature for indels).
5. Build a phylogenetic clone tree based on all filtered mutations.

This pipeline is still in active development. The full set of features that we *hope* to implement is illustrated below.

![workflow](docs/images/workflow.png)

## Dependencies

- [Nextflow](https://www.nextflow.io/)
- [hairpin](https://github.com/cancerit/hairpin-wrapper) (wrapper by CASM IT at the Wellcome Sanger institute)
- [vcfilter]() (wrapper by CASM IT at the Wellcome Sanger institute, public repo coming soon)
- [tabix](https://www.htslib.org/doc/tabix.html)
- [cgpVAFcommand](https://github.com/cancerit/vafCorrect/tree/dev)
- [R](https://www.r-project.org/)
- [MutationsPy](https://github.com/Phuong-Le/MutationsPy) (no need to install independently as it is integrated as a git submodule)

## Installation

Clone this repository, including the MutationsPy submodule

```
git clone --recursive git@github.com:Phuong-Le/process_sanger_lcm-nf.git
```

## Usage

### samplesheet.csv

First, prepare a samplesheet with your input data. Demo samplesheets for SNVs and indels can be found in the `demo_files/` directory. The samplesheet must be in CSV format and should look as follow:

| sample_id       | match_normal_id | pdid    | bam_match                    | bai_match                        | vcf                                                 | vcf_tbi                                                 | bam                                 | bai                                     |
| --------------- | --------------- | ------- | ---------------------------- | -------------------------------- | --------------------------------------------------- | ------------------------------------------------------- | ----------------------------------- | --------------------------------------- |
| PD47151n_lo0002 | PD47151b        | PD47151 | path/to/PD47151/PD47151b.bam | path/to/PD47151/PD47151b.bam.bai | path/to/PD47151/PD47151n_lo0002.pindel.annot.vcf.gz | path/to/PD47151/PD47151n_lo0002.pindel.annot.vcf.gz.tbi | path/to/PD47151/PD47151n_lo0002.bam | path/to/PD47151/PD47151n_lo0002.bam.bai |
| PD47151n_lo0004 | PD47151b        | PD47151 | path/to/PD47151/PD47151b.bam | path/to/PD47151/PD47151b.bam.bai | path/to/PD47151/PD47151n_lo0004.pindel.annot.vcf.gz | path/to/PD47151/PD47151n_lo0004.pindel.annot.vcf.gz.tbi | path/to/PD47151/PD47151n_lo0004.bam | path/to/PD47151/PD47151n_lo0004.bam.bai |
| PD52103n_lo0002 | PD52103b        | PD52103 | path/to/PD52103/PD52103b.bam | path/to/PD52103/PD52103b.bam.bai | path/to/PD52103/PD52103n_lo0002.pindel.annot.vcf.gz | path/to/PD52103/PD52103n_lo0002.pindel.annot.vcf.gz.tbi | path/to/PD52103/PD52103n_lo0002.bam | path/to/PD52103/PD52103n_lo0002.bam.bai |

Each row represents one sample. Sample IDs should be unique.

N.B. VCF files must gzipped, with the extension `*.vcf.gz`. Additionally, due to `cgpVAF` requirements, the BAM files (`bam`, `bam_match`) should be in the form `<sample_id>.bam `and BAI files (`bai`, `bai_match`) in the form `<sample_id>.bam.bai`. These constraints may be lifted soon (see **Future work**).

### Running the pipeline

Now, you can run the pipeline using:

```
nextflow run process_sanger_lcm-nf \
    --samplesheet /path/to/samplesheet.csv \
    --mut_type snv \
    --reference_genome /path/to/genome.fa \
    --high_depth_bed /path/to/HiDepth.bed.gz \
    --outdir out/
```

### Running the pipeline on Sanger's Farm

If you are running this pipeline from the Wellcome Sanger Institute HPC cluster ('the farm'), you can run the pipeline as follows:

```
# load modules
module load nextflow/23.10.1-5891
module load hairpin/1.0.6
module load vcfilter/1.0.4
module load tabix/1.13
module load cgpVAFcommand/2.5.0
module load R/4.1.0

# run the pipeline from a head or interactive node
nextflow run /path/to/process_sanger_lcm-nf \
  -c conf/sanger_lsf.config \
  --samplesheet samplesheet.csv \
  --mut_type snv \
  --reference_genome /path/to/reference_data/hg38/genome.fa \
  --high_depth_bed /path/to/reference_data/hg38/highdepth.bed.gz \
  --outdir out/

# run the pipeline from a compute node
bsub \
  -cwd "./" \
  -o log/%J.${mut_type}.out \
  -e log/%J.${mut_type}.err \
  -R "select[mem>5000] rusage[mem=5000]" -M5000 \
  "nextflow run /path/to/process_sanger_lcm-nf -c conf/sanger_lsf.config --samplesheet samplesheet.csv --mut_type snv --reference_genome /path/to/reference_data/hg38/genome.fa --high_depth_bed /path/to/reference_data/hg38/highdepth.bed.gz --outdir out/"

```

N.B. This pipeline (particularly the sanger_lsf.config) has only be tested on non-malignant colon data (~2,000-9,000 mutations per sample). If you have a bigger dataset, it might be neccessary to experiment with different memory settings.

### All parameters

In order to view all available parameters, their types, descriptions and defaults, please use the `--help` argument.

```
$ nextflow run process_sanger_lcm-nf --help 

 N E X T F L O W   ~  version 24.04.4

Launching `./main.nf` [boring_lumiere] DSL2 - revision: 91dc7941f3

Typical pipeline command:

  nextflow run process_sanger_lcm-nf \
    --samplesheet /path/to/samplesheet.csv \
    --mut_type snv \
    --reference_genome /path/to/genome.fa \
    --high_depth_bed /path/to/HiDepth.bed.gz \
    --outdir out/

Input/output options
  --samplesheet                     [string]  A file that contains paths to relevant files for the processes in the pipeline, with headers similar to 
                                               demo_files/lcm_processing_input_indel.tsv for Indels or demo_files/lcm_processing_input_snv.tsv for SNVs. 
  --mut_type                         [string]  The type of mutation in the input files. (accepted: snv, indel) [default: snv]
  --reference_genome                 [string]  Path to the reference genome fasta file.
  --with_match_normal                [boolean] null [default: true]
  --outdir                           [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud 
                                               infrastructure. 
  --help                             [boolean] Return pipeline options and example command.

Conpair options
  --conpair                          [boolean] Run Conpair? [default: true]
  --marker_bed                       [string]  null
  --marker_txt                       [string]  null
  --concordance_threshold            [integer] null [default: 90]
  --contamination_threshold_samples  [number]  null [default: 0.3]
  --contamination_threshold_match    [integer] null [default: 5]

SNV filtering options
  --filter_snv                       [boolean] Filter SNVs? [default: true]

Indel filtering options
  --filter_indel                     [boolean] Filter indels? [default: true]

vcfilter options
  --vcfilter_config                  [string]  Pre-cgpVAF filtering. Defaults to data/indel_default.filter for Indels and data/snv_default.filter for 
                                               SNVs. 

cgpVAF options
  --high_depth_bed                   [string]  Path to a reference BED file of high depth regions.

Phylogenetics options
  --phylogenetics                    [boolean] Build a phylogenetic tree from the mutations? [default: true]
  --with_topology                    [string]  null
  --phylogenetics_outdir_basename    [string]  null

```

## Authors

[Phuong Le](https://github.com/Phuong-Le) (al35@sanger.ac.uk)

[Rashesh Sanghvi](https://github.com/Rashesh7) (rs30@sanger.ac.uk)

[Alex Tidd](https://github.com/alextidd) (at31@sanger.ac.uk)

Chuling Ding (cd23@sanger.ac.uk)

## Future work

* [ ] Improve cgpVAF runtime by parallelising computation by chromosome and/or batches of samples.
* [ ] Resolve the constraints on the file names in the input (see N.B. in `samplesheet.csv` section).
* [ ] Update hairpin to hairpin2.
* [ ] Containerise all steps so that the pipeline will run off-farm.
* [ ] Add VCF reflagging.
* [ ] Change from cgpVAF to bamtoR-based VAF calculation.

## Acknowledgements

Thanks to Chloe Pacyna and Yichen Wang for guiding me through the necessary steps to process LCM data.

Thanks to Shriram Bhosle from CASM IT for showing me how to use cgpVaf (ie vafCorrect).

The BetaBinomial filtering of germline and artefacts was adapted from [code written by Tim Coorens](https://github.com/TimCoorens/Unmatched_NormSeq/tree/master).
