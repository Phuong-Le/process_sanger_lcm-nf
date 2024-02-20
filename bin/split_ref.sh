#!/usr/bin/env bash


# parsing arguments
while [[ "$#" -gt 0 ]]
do
  case "$1" in
    --reference )
      reference="$2"
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



seq=$(grep '>' ${reference} | awk '{print $1}' | sed 's/^>//')
mv $reference "${seq}.fa"