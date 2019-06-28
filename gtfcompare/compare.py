#!/usr/bin/env python
# -*- coding: utf-8 -*-

from argument_parser import argument_parser
from to_gtftree import to_gtftree
import time

reference, input_files = argument_parser()

time_begin = time.time()
ref_chrom_dict, ref_introns = to_gtftree(reference, 'ref')
time_end = time.time()
print("ref", time_end - time_begin)
for input_file in input_files:
    time_begin = time.time()
    quer_chrom_dict, quer_introns = to_gtftree(input_file, input_file.split('/')[-1])
    time_end = time.time()
    print(input_file.split('/')[-1], time_end - time_begin)