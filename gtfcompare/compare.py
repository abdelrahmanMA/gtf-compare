#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
from argument_parser import argument_parser
from to_gtftree import to_gtftree
from interval_matching import match_intervals
from writers import write_exons, write_itrack

time_very_begin = time.time()

reference, input_files = argument_parser()

time_begin = time.time()
ref_chrom_dict, ref_introns = to_gtftree(reference, 'ref')
time_end = time.time()
print("ref", time_end - time_begin)
for input_file in input_files:
    shortname = input_file.split('/')[-1]
    time_begin = time.time()
    quer_chrom_dict, quer_introns = to_gtftree(input_file, shortname)
    time_end = time.time()
    print(shortname, time_end - time_begin)
    time_begin = time.time()
    match_intervals(quer_chrom_dict, ref_chrom_dict)
    time_end = time.time()
    print('match intervals', time_end - time_begin)
    time_begin = time.time()
    write_exons(shortname, quer_chrom_dict)
    time_end = time.time()
    print('write dotexon', time_end - time_begin)
    time_begin = time.time()
    write_exons('ref_' + shortname, ref_chrom_dict)
    time_end = time.time()
    print('write dotexon ref', time_end - time_begin)
    time_begin = time.time()
    write_itrack(shortname, quer_chrom_dict, ref_introns)
    time_end = time.time()
    print('write itrack', time_end - time_begin)
    time_begin = time.time()
    write_itrack('ref_' + shortname, ref_chrom_dict, quer_introns)
    time_end = time.time()
    print('write itrack ref', time_end - time_begin)
print('finished in', time.time() - time_very_begin)