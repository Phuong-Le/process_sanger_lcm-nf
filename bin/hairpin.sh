#!/usr/bin/env bash

# parsing arguments
while [[ "$#" -gt 0 ]]
do
  case "$1" in
    --vcf )
      vcf="$2"
      shift 2
      ;;
    --bam )
      bam="$2"
      shift 2
      ;;
    --vcfilter_config )
      vcfilter_config="$2"
      shift 2
      ;;
    --)
      shift;
      break
      ;;
    *)
      echo "Unexpected option: $1"
      ;;
  esac
done


# hairpin annotation
hairpin \
    -v $vcf \
    -b $bam \
    -o . \
    -g 38
vcf_hairpin=$(basename $vcf | sed 's/vcf.gz/hairpin.vcf/') # output file from hairpin 

# vcfilter
vcfilter filter -o . -i ${vcfilter_config} ${vcf_hairpin}
vcf_hairpin_filter=$(basename $vcf | sed 's/vcf.gz/hairpin.filter.vcf/') # output file from filter

# prepare for cgpVaf, might move to cgpVaf
bgzip $vcf_hairpin_filter
vcf_hairpin_filter_bgzipped=$(basename $vcf | sed 's/vcf.gz/hairpin.filter.vcf.gz/')
tabix $vcf_hairpin_filter_bgzipped 

#