#!/usr/local/bin/python

##   LICENSE
# Copyright (c) 2018-2019 Genome Research Ltd.
# Author: Cancer Genome Project cgpit@sanger.ac.uk
#
#
# This file is part of vcf_flag_modifier.
#
# CaVEMan is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#    1. The usage of a range of years within a copyright statement contained within
#    this distribution should be interpreted as being equivalent to a list of years
#    including the first and last year specified and all consecutive years between
#    them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
#    2009, 2011-2012’ should be interpreted as being identical to a statement that
#    reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
#    statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
#    identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
#    2009, 2010, 2011, 2012’."
#
#

"""
Main script for the flag_modifier project.
A python tool to modify flags in a VCF file.
Offers options of a list of flags, or restricting to a bed file of
locations, even to a bed file of locations with specific changes
"""
# Core imports
import pkg_resources  # part of setuptools
import argparse
import os
import sys

from vcfflagmodifier.VcfParse import VcfParse

VERSION = pkg_resources.require("vcfflagmodifier")[0].version


def generate_arg_str(parsed_args):
    ag_str = ""
    idx = 0
    for arg in sorted(parsed_args.__dict__.keys()):
        val = parsed_args.__dict__[arg]
        if not val:
            continue
        if arg == "vcfin" or arg == "vcfout" or arg == "bedfile":
            val = os.path.basename(val)

        if arg == "flagremove":
            inner_idx = idx
            for v in val:
                if inner_idx > 0:
                    ag_str += ","
                ag_str += "{}={}".format(arg, v)
                inner_idx += 1
        else:
            if idx > 0:
                ag_str += ","
            ag_str += "{}={}".format(arg, val)
        idx += 1
    return ag_str


parser = argparse.ArgumentParser(
    description="Modify a VCF file, removing provided flags.",
    prog="vcf_flag_modifier.py",
)

parser.add_argument(
    "-f",
    "--vcfin",
    dest="vcfin",
    metavar="input.vcf[.gz]",
    help="Path to input VCF file",
    required=True,
)

parser.add_argument(
    "-o",
    "--vcfout",
    metavar="output.vcf",
    help="path to write output VCF file",
    required=False,
    default="out.vcf",
)

parser.add_argument(
    "-b",
    "--bedfile",
    metavar="positions.bed[.gz]",
    dest="bedfile",
    help="""Path to bed file of positions to consider for flag modification
                        Bed file can be `contig  start  stop` or can include specific changes for each coordinate
                        `contig  start  stop    ref alt`""",
    required=False,
)

parser.add_argument(
    "-a",
    "--alleles",
    dest="alleles",
    action="store_true",
    help="""Path to bed file of positions to consider for flag modification
                        Bed file can be `contig  start  stop` or can include specific changes for each coordinate
                        `contig  start  stop    ref alt`""",
    required=False,
)

parser.add_argument(
    "-v",
    "--version",
    dest="version",
    action="version",
    help="Print version information",
    version="%(prog)s " + VERSION,
)

list_grp = parser.add_mutually_exclusive_group()

list_grp.add_argument(
    "-t",
    "--listflags",
    dest="listflags",
    action="store_true",
    help="Instead of modifying a VCF, list flags avaialble from header",
    required=False,
    default=False,
)
list_grp.add_argument(
    "-l",
    "--flagremove",
    dest="flagremove",
    metavar="UMV",
    action="append",
    help="Flag(s) to remove from the passed VCF file (repeatable)",
    required=False,
)

args = parser.parse_args()

arg_str = generate_arg_str(args)
vcfparse = None

if args.bedfile:
    vcfparse = VcfParse(
        args.vcfin,
        args.vcfout,
        os.path.basename(__file__),
        arg_str,
        bedfile=args.bedfile,
        alleles=args.alleles,
    )
else:
    vcfparse = VcfParse(args.vcfin, args.vcfout, os.path.basename(__file__), arg_str)

vcfparse.get_header_flags()

if args.listflags is True:
    vcfparse.print_header_flags()
else:
    vcfparse.check_filters_against_flaglist(args.flagremove)
    vcfparse.refilter_vcf()
    print("Refiltering complete")
