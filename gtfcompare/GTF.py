#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
class Feature():
    '''
    A Class to represent a GTF feature/line
    '''
    def __init__ (self, chromosome, source, f_type, start, end, score, strand, phase, attributes):
        self.chromosome = chromosome
        self.chrom = self.chromosome
        self.source = source
        self.f_type = f_type
        self.start = start
        self.begin = self.start
        self.end = end
        self.size = end - start + 1
        self.score = score
        self.strand = strand
        self.phase = phase
        self.frame = self.phase
        self.attributes = attributes
        self.gene_id = attributes.get('gene_id', None)
        self.gene_name = attributes.get('gene_name', None)
        self.transcript_id = attributes.get('transcript_id', None)

    def set_id(self):
        if self.f_type == 'exon':
            self.id = self.attributes['exon_id']
        elif self.f_type == 'transcript' or self.f_type == 'mRNA':
            self.id = self.attributes['transcript_id']
        elif self.f_type == 'gene':
            self.id = self.attributes['gene_id']

class GTF():
    '''
    Class to repreent GTF file
    '''
    def __init__(self, GTF) :
        self.GTF = GTF

    def order(self, exons, exon_id):
        re_len = len(exons)
        for i, ex in enumerate(exons):
            ex.attributes['exon_number'] = re_len - i
            ex.attributes['exon_id'] = exon_id - i - 1
            ex.set_id()
            yield ex

    def parse(self):
        # Open the file with read only permit
        GTF = open(self.GTF)
        # use readline() to read the first line
        line = GTF.readline()
        # use the read line to read further.
        # If the file is not empty keep reading one line
        # at a time, till the file is empty
        exon_id = 1
        exon_number = 1
        neg_ordered = False
        re_order = []
        prev_transcript = None
        prev_start = None
        gff3 = False
        gff3_map = {}
        while line:
            if line[0:1] != '#':
                # Read each line into 9 columns
                columns = line.split('\t')
                chromosome = columns[0]
                source = columns[1]
                f_type = columns[2]
                start = int(columns[3])
                end = int(columns[4])
                score = columns[5]
                strand = columns[6]
                phase = columns[7]
                attributes = columns[8]
                if gff3:
                    r1 = re.findall(r'(ID=[^(;|"|\n)]*|geneID=[^(;|"|\n)]*|gene_name=[^(;|"|\n)]*)', attributes)
                    if r1:
                        if f_type == 'mRNA':
                            gff3_map[r1[0].split('=')[1]] = {'gene_id':r1[1].split('=')[1], 'transcript_id': r1[0].split('=')[1], 'gene_name':r1[2].split('=')[1]}
                    r1 = re.findall(r'(Parent[^(;|"\n)]*)', attributes)
                    if r1:
                        needed_attributes = gff3_map[r1[0].split('=')[1]]
                else:
                    r1 = re.findall(r'(gene_id "[^(;|"|\n)]*|transcript_id "[^(;|"|\n)]*|gene_name "[^(;|"|\n)]*)', attributes)
                    needed_attributes = {x.split(' "')[0]:x.split(' "')[1] for x in r1}
                if f_type == 'exon':
                    exon = Feature(chromosome, source, f_type, start, end, score, strand, phase, needed_attributes)
                    current_transcript = needed_attributes['transcript_id']
                    if strand == '+' or strand == '.' or neg_ordered:
                        if neg_ordered and re_order:
                            for i, re_ord in enumerate(re_order):
                                re_ord.attributes['exon_number'] = i + 1
                                exon_number = i + 1
                                re_ord.attributes['exon_id'] = exon_id - len(re_order) + i
                                re_ord.set_id()
                                yield re_ord
                            re_order = []
                        if prev_transcript == current_transcript:
                            exon_number += 1
                        else:
                            exon_number = 1
                        exon.attributes['exon_number'] = exon_number
                        exon.attributes['exon_id'] = exon_id
                        exon.set_id()
                        yield exon

                    else:
                        if prev_transcript == current_transcript or prev_transcript is None:
                            if prev_start is None:
                                re_order.append(exon)
                            elif prev_start < end:
                                re_order.append(exon)
                            else:
                                exon.attributes['exon_number'] = 1
                                exon.attributes['exon_id'] = exon_id
                                re_order.append(exon)
                                neg_ordered = True
                        else:
                            for ex in self.order(re_order, exon_id):
                                yield ex
                            re_order = [exon]
                    exon_id += 1
                    prev_transcript = current_transcript
                    prev_start = start
                elif re_order:
                    for ex in self.order(re_order, exon_id):
                        yield ex
                    re_order = []
            else:
                line = line[2:]
                if 'gff-version' in line:
                    gffv = line.replace('gff-version ', '').rstrip()
                    if  gffv == '3':
                        gff3 = True
            # use realine() to read next line
            line = GTF.readline()
        if re_order:
            for ex in self.order(re_order, exon_id):
                yield ex
        GTF.close()