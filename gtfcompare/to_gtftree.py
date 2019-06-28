#!/usr/bin/env python
# -*- coding: utf-8 -*-

from intervaltree import Interval, IntervalTree
import numpy as np
from gtftree import GtfInterval, GtfExon, GtfTranscript, GtfGene
from GTF import GTF

def to_gtftree(Gtf_file, shortname):
    """Generates GTF stats and constructs GTF dictionaries"""
    GTF_exons = GTF(Gtf_file)
    gtf_chrom_dict = {}

    gtf_introns = IntervalTree()  # Interval Tree to hold reference introns
    # Next line defines a variable which will be used to calculate exonic unique length
    gtf_ulen = {'+': IntervalTree(), '-': IntervalTree(), '.': IntervalTree()}
    # Next line defines a variable which will be used to calculate reference statistics
    GTF_STATS = {'Total_Exons': [], 'Unique_Exons': [], '5_Prime_UTR': 0, '3_Prime_UTR': 0, 'Intronic_Edges': 0}

    gene_strand_count = {}
    tran_strand_count = {}
    gtf_chrom_dict = {}
    gtf_gene_dict = {}
    gtf_transcript_dict = {}
    gene_name_id = {}
    number = 0
    intron_start = 0
    intron_end = 0
    last_interval = None
    prev_transcript_id = None
    prev_strand = None
    reverse = None

    # Loops over all reference exons
    for i, exon in enumerate(GTF_exons.parse()):
        exon_interval = Interval(exon.begin, exon.end + 1)

        # Calculates Unique length
        mn = exon.begin
        mx = exon.end + 1
        # Loops over all overlapping intervals removes them
        for b, e, d in gtf_ulen[exon.strand][mn - 1:mx]:
            if b < mn:
                mn = b
            if e > mx:
                mx = e
            gtf_ulen[exon.strand].removei(b, e, d)
        # Adds and interval to unique length starting from smallest start to the largest end
        gtf_ulen[exon.strand].addi(mn, mx, i)

        # Checks if the current exon in a new transcript
        if prev_transcript_id != exon.transcript_id:
            if last_interval:
                if number == 'First':
                    last_interval.transcriptIds[prev_transcript_id] = 'Single'
                    GTF_STATS['3_Prime_UTR'] += 1
                elif prev_strand == '-' and reverse:
                    last_interval.transcriptIds[prev_transcript_id] = 'First'
                    GTF_STATS['5_Prime_UTR'] += 1
                else:
                    last_interval.transcriptIds[prev_transcript_id] = 'Last'
                    GTF_STATS['3_Prime_UTR'] += 1

                GTF_STATS['Intronic_Edges'] -= 1

            if exon.attributes['exon_number'] == 1:
                number = 'First'
            else:
                number = 'Last'

            if exon.strand == '-' and number == 'Last':
                reverse = True

            prev_transcript_id = exon.transcript_id
            prev_strand = exon.strand

            if (exon.strand == '-' and not reverse):
                intron_start = None
                intron_end = exon.start
            else:
                intron_start = exon.end
                intron_end = None
        else:
            number = exon.attributes['exon_number']

            if (exon.strand == '-' and not reverse):
                intron_start = exon.end
            else:
                intron_end = exon.start

        GTF_STATS['Total_Exons'].append(exon.size)
        if intron_start and intron_end:
            gtf_introns[intron_start:intron_end] = exon.strand

        if (exon.strand == '-' and not reverse):
            intron_end = exon.start
        else:
            intron_start = exon.end

        if number == 'Single':
            exit(1)
        elif number == 'First':
            GTF_STATS['5_Prime_UTR'] += 1
            GTF_STATS['Intronic_Edges'] += 1
        elif number == 'Last':
            GTF_STATS['Intronic_Edges'] += 1
            GTF_STATS['3_Prime_UTR'] += 1
        else:
            GTF_STATS['Intronic_Edges'] += 1
            GTF_STATS['Intronic_Edges'] += 1

        # To Replace gene id with gene name
        if exon.gene_name:
            gene_name_id[exon.gene_id] = exon.gene_name

        # Creates a dictionary for the chromosome if it didn't exist
        if exon.chrom not in gtf_chrom_dict:
            gtf_chrom_dict[exon.chrom] = {}

        # Creates a dictionary for the strand if it didn't exist
        if exon.strand not in gtf_chrom_dict[exon.chrom]:
            # Creates an interval tree for intervals and exons dictionary relative to the current chromosome and strand
            gtf_intervals = IntervalTree()
            gtf_exons_dict = {}

            last_interval = GtfInterval(exon_interval, exon, number)
            gtf_intervals[exon.start: exon.end + 1] = last_interval

            last_exon = GtfExon(exon.id, last_interval, exon.transcript_id, exon.gene_id)
            gtf_exons_dict[exon.id] = last_exon

        else:
            # Retrives the interval tree and exon dictionary
            gtf_intervals, gtf_exons_dict = gtf_chrom_dict[exon.chrom][exon.strand]
            found = False
            dz = None
            # Checks if the interval already exists
            for b, e, z in gtf_intervals.envelop(exon.begin, exon.end + 1):
                if b == exon_interval.begin and e == exon_interval.end:
                    found = True
                    dz = z
                    break
            # Means the interval didn't exist
            if not found:
                last_interval = GtfInterval(exon_interval, exon, number)
                gtf_intervals[exon.begin: exon.end + 1] = last_interval
                # GTF_STATS['Unique_Exons'].append(exon.size)
            else:
                if dz is None:
                    exit(2)
                # Add exon data to the existing interval
                dz.add(exon.id, exon.transcript_id, number, exon.gene_id)
                last_interval = dz

            if exon.id not in gtf_exons_dict:
                last_exon = GtfExon(exon.id, last_interval, exon.transcript_id, exon.gene_id)
                gtf_exons_dict[exon.id] = last_exon
            else:
                print('Duplicate Exon ID', exon.id, exon_interval)
                exit(3)

        gtf_chrom_dict[exon.chrom][exon.strand] = [gtf_intervals, gtf_exons_dict]

        if exon.transcript_id not in gtf_transcript_dict:
            last_transcript = GtfTranscript(exon.transcript_id, last_exon, exon.gene_id)
            gtf_transcript_dict[exon.transcript_id] = last_transcript
        else:
            last_transcript = gtf_transcript_dict[exon.transcript_id]
            if last_transcript.gene_id != exon.gene_id:
                print('Duplicate Transcript id in different genes', last_transcript.id)
                exit(4)
            elif last_transcript.id != exon.transcript_id:
                print('Duplicate Transcript', last_transcript.id)
                exit(4)
            last_transcript.add(last_exon)

        if exon.gene_id not in gtf_gene_dict:
            gtf_gene_dict[exon.gene_id] = GtfGene(exon.gene_id, last_transcript)
        else:
            if gtf_gene_dict[exon.gene_id].id != exon.gene_id:
                print('Duplicate Gene', gtf_gene_dict[exon.gene_id].id)
                exit(5)
            gtf_gene_dict[exon.gene_id].add(last_transcript)

        # Next few line sets the strand count for each gene and transcript to determine their strand later
        if exon.chrom not in tran_strand_count:
            tran_strand_count[exon.chrom] = {}
        if exon.transcript_id not in tran_strand_count[exon.chrom]:
            tran_strand_count[exon.chrom][exon.transcript_id] = {'+': 0, '-': 0, '.': 0}
        tran_strand_count[exon.chrom][exon.transcript_id][exon.strand] += 1

        if exon.chrom not in gene_strand_count:
            gene_strand_count[exon.chrom] = {}
        if exon.gene_id not in gene_strand_count[exon.chrom]:
            gene_strand_count[exon.chrom][exon.gene_id] = {'+': 0, '-': 0, '.': 0}
        gene_strand_count[exon.chrom][exon.gene_id][exon.strand] += 1

    if number == 'First':
        last_interval.transcriptIds[prev_transcript_id] = 'Single'
        GTF_STATS['3_Prime_UTR'] += 1
    elif exon.strand == '-' and reverse:
        last_interval.transcriptIds[prev_transcript_id] = 'First'
        GTF_STATS['5_Prime_UTR'] += 1
    else:
        last_interval.transcriptIds[prev_transcript_id] = 'Last'
        GTF_STATS['3_Prime_UTR'] += 1

    GTF_STATS['Intronic_Edges'] -= 1

    for chromo in gtf_chrom_dict:
        for strand in gtf_chrom_dict[chromo]:
            gtf_chrom_dict[chromo][strand].append({})
            gtf_chrom_dict[chromo][strand].append({})
    # Determines each transcript's strand
    for chromo in tran_strand_count:
        for tran, strands in tran_strand_count[chromo].items():
            if strands['+'] > strands['.'] and strands['+'] > strands['-']:
                strand = '+'
            elif strands['-'] > strands['.'] and strands['-'] > strands['+']:
                strand = '-'
            else:
                strand = '.'
            gtf_chrom_dict[chromo][strand][2][tran] = gtf_transcript_dict[tran]
    # Determines each gene's strand
    for chromo in gene_strand_count:
        for gene, strands in gene_strand_count[chromo].items():
            if strands['+'] > strands['.'] and strands['+'] > strands['-']:
                strand = '+'
            elif strands['-'] > strands['.'] and strands['-'] > strands['+']:
                strand = '-'
            else:
                strand = '.'
            gtf_chrom_dict[chromo][strand][3][gene] = gtf_gene_dict[gene]

    # Write temporarily Reference stats file
    with open(shortname + '.tstat', 'w') as gtfstatf:
        wf = False
        for st in ['+', '-', '.']:
            for b, e, d in gtf_ulen[st]:
                GTF_STATS['Unique_Exons'].append(e - b)

        for key in ['Total_Exons', 'Unique_Exons']:
            value = GTF_STATS[key]
            val = np.array(value)
            gtfstatf.write('{}(count|total|mean|std)\t{}|{}|{}|{}\n'.
                        format(key, len(val), val.sum(), round(val.mean(), 1), round(val.std(), 1)))

        for key in ['5_Prime_UTR', '3_Prime_UTR', 'Intronic_Edges']:
            value = GTF_STATS[key]
            if key == '5_Prime_UTR':
                if value == GTF_STATS['3_Prime_UTR']:
                    continue
                else:
                    print('Something went wrong with stats UTRS don\'t match')
                    wf = True
            if key == '3_Prime_UTR' and not wf:
                key = 'Transcripts'
            gtfstatf.write('{}(count)\t{}\n'.format(key, value))
    return [gtf_chrom_dict, gtf_introns]
