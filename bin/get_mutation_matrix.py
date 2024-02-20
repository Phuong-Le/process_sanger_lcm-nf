#!/usr/bin/env python3

from collections import Counter
import pandas as pd

from mutationsPy.read_file import read_vcf, extract_seq_from_fasta
from mutationsPy.get_mut import get_context
from mutationsPy.gen_context import gen_context

import argparse, sys

def get_mutation_matrix(ref_dir, vcf_path, sample_id, outfile, kmer=3, sep="\t"):
    vcf = read_vcf(vcf_path, usecols = ['CHROM', 'POS', 'REF', 'ALT'])
    chroms = set(vcf['CHROM'])
    for chrom in chroms:
        ref_path = f'{ref_dir}/{chrom}.fa'
        seq = extract_seq_from_fasta(ref_path, identifier=chrom)
        mut_count = Counter(symmetrise_context(get_context(seq, pos, ref, alt)) for pos, ref, alt in zip(vcf_chrom['POS'], vcf_chrom['REF'], vcf_chrom['ALT']))
        mut_counts.update(mut_count)
    mutations = gen_context(kmer)
    mut_counts_vector = [str(mut_counts[m]) if m in mut_counts else 0 for m in mutations]
    
    # write file to disk
    header = "sample_id" + sep + sep.join(mutations) 
    content = sample_id + sep + sep.join(mut_counts_vector)
    with open(outfile, mode='wt') as file:
        file.writelines(f'{header}\n{content}\n')
 
def get_arguments():
    parser = argparse.ArgumentParser(description='VCF to mutation matrix')
    parser.add_argument('--ref_dir', required=True,
                        help='directory that contains the reference sequence(s)', type = str)
    parser.add_argument('--vcf_path', required=True,
                        help='path to vcf file to be converted', type = str)
    parser.add_argument('--sample_id', required=True,
                        help='sample ID', type = str)
    parser.add_argument('--kmer', '-k', required=False,
                        help='kmer size, default to 3', type = int)
    parser.add_argument('--sep', required=False,
                        help='delimiter of output file, default to tab', type = str)
    parser.add_argument('--outfile', '-o', required=True,
                        help='output file', type = str)
    return parser

def main(args):
    get_mutation_matrix(arg.ref_dir, arg.vcf_path, arg.sample_id, arg.outfile, kmer=arg.kmer, sep=arg.sep)

if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))