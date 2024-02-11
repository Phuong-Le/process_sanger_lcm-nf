#!/usr/bin/env bash

# parsing arguments
while [[ "$#" -gt 0 ]]
do
  case "$1" in
    --vcf )
      vcf="$2"
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



# vcfilter
vcfilter filter -o . -i ${vcfilter_config} ${vcf}
vcf_filtered=$(basename $vcf | sed 's/vcf.gz/filter.vcf/') # output file from filter

# prepare for cgpVaf, might move to cgpVaf
bgzip $vcf_filtered
vcf_filtered_bgzipped=$(basename $vcf | sed 's/vcf.gz/filter.vcf.gz/')
tabix $vcf_filtered_bgzipped 