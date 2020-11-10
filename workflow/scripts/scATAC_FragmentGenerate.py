# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-01-17 01:43:08
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-06-10 23:32:28

import argparse
import pysam
import time, os
import sys
from collections import defaultdict

chr_list = ['chr'+str(i) for i in list(range(1,23))] 
chr_list = chr_list + ['chrX', 'chrY']

def CommandLineParser():
    parser=argparse.ArgumentParser(description = "This is a description of input args")
    parser.add_argument("-b","--bam", dest = "bamfile", default = "",help = "The merged bamfile.")
    # parser.add_argument("-b","--barcodefq",dest = "barcode_fastq",default = "",help = "The barcode fastq file (for sciATAC or 10x genomics, R2)")
    parser.add_argument("-o", "--outfile", nargs = "?", default=sys.stdout, type = argparse.FileType('w'), help = "The output file.")
    parser.add_argument("--addtag", dest = "addtag", default = "", help = "The tag cell barcodes will be added to. If not set, no tag will be added.")
    parser.add_argument("--cell_barcode", required=True, help = "Where the cell identities are stored. For example, 'CB'. "
        "If not set, cell identities will be inferred from read name.")
    return parser.parse_args()

bamfile_in = pysam.AlignmentFile(parser.bamfile,"rb")

for read in bamfile_in:
    if read.has_tag(parser.cell_barcode):
        barcode = read.get_tag(parser.cell_barcode)
        flag = read.flag
        if read.reference_name in chr_list and (flag & 0x2 !=0) and (flag & 0xc == 0) and (flag & 0x900 == 0) and read.mapping_quality > 30 and read.template_length < 1000 and read.template_length > 10:
            frag_list = [read.reference_name, str(read.reference_start + 4), str(read.reference_start + read.template_length - 5), barcode]
            print(*frag_list, sep = '\t', file = parser.outfile)