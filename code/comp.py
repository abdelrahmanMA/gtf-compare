import gffutils
from argparse import ArgumentParser
from intervaltree import Interval, IntervalTree
import time
from tqdm import tqdm
import itertools
import numpy as np
import six


time_begin = time.time()
input_files = None
gene_name_id = {}  # gene name dictionary
ref_introns = IntervalTree()  # Interval Tree to hold reference introns

interval_matches = {}  # dictionary that contains all query matches
interval_best_matches = {}  # dictionary that contains query best matches
ref_matches = {}  # dictionary that contains all reference matches
ref_best_matches = {}  # dictionary that contains reference best matches

NOTES = ['-',
         'A match was found on the unknown strand of the reference',
         'A match was found on the opposite strand on the reference',
         'A match was found on the + and - strands of the reference query\'s strand is unknown',
         'More than one best match was found all on the same strand same as query strand',
         'More than one best match was found on the unknown strand of reference',
         'More than one best match was found all on the same strand opposite of query strand',
         'More than one best match was found on both strands of reference query\'s strand is unknown',
         '']


class Exon(gffutils.Feature):
    """A Class To extend gffutils to add size feature and to make it easier to access gene id and transcript id"""

    def __init__(self, feature):
        gffutils.Feature.__init__(self, feature.seqid, feature.source, feature.featuretype,
                                  feature.start, feature.end, feature.score, feature.strand, feature.frame,
                                  feature.attributes, feature.extra, feature.bin, feature.id, feature.dialect,
                                  feature.file_order, feature.keep_order, feature.sort_attribute_values)
        self.size = feature.end - feature.start + 1
        # if feature.attributes.get('gene_name'):
        #     self.gene_id = feature.attributes.get('gene_name')[0]
        # else:
        self.gene_id = feature.attributes.get('gene_id')[0]
        self.transcript_id = feature.attributes.get('transcript_id')[0]
        self.begin = feature.start


class GtfInterval:
    """Represents a GTF Interval"""

    def __init__(self, interval, exon, number):

        self.interval = interval
        self.begin = self.interval.begin
        self.end = self.interval.end
        self.exonIds = {exon.id}
        self.transcriptIds = {exon.transcript_id: number}
        self.geneIds = {exon.gene_id}
        self.strand = exon.strand
        self.chrom = exon.chrom
        self.note = None

    def add(self, exon_id, transcript_id, number, gene_id):
        """Adds all possible exon ids and transcript ids and gene ids"""

        if exon_id not in self.exonIds:
            self.exonIds.add(exon_id)
        else:
            print('An exon id was found twice')

        if transcript_id not in self.transcriptIds:
            self.transcriptIds[transcript_id] = number
        else:
            print('An transcript was found twice')

        if gene_id not in self.geneIds:
            self.geneIds.add(gene_id)


class GtfExon:
    """Represents GTF Exon"""

    def __init__(self, _id, gtf_interval, transcript_id, gene_id):

        self.id = _id
        self.gtf_interval = gtf_interval
        self.begin = self.gtf_interval.begin
        self.end = self.gtf_interval.end
        self.strand = self.gtf_interval.strand
        self.chrom = self.gtf_interval.chrom
        self.transcript_id = transcript_id
        self.gene_id = gene_id


class GtfTranscript:
    """Represents GTF Transcript"""

    def __init__(self, _id, gtf_exon, gene_id):
        interval = gtf_exon.gtf_interval.interval
        self.id = _id
        self.gtf_exons = {gtf_exon}
        self.exonIds = {gtf_exon.id}
        self.gene_id = gene_id
        self.begin = interval.begin
        self.end = interval.end

    def add(self, gtf_exon):
        """Adds Exons to GTF Transcript"""
        interval = gtf_exon.gtf_interval.interval
        self.gtf_exons.add(gtf_exon)
        self.exonIds.add(gtf_exon.id)
        self.begin = min(self.begin, interval.begin)
        self.end = max(self.end, interval.end)


class GtfGene:
    """Represents GTF Gene"""

    def __init__(self, _id, gtf_transcript):
        self.id = _id
        self.gtf_transcripts = {gtf_transcript}
        self.transcriptIds = {gtf_transcript.id}
        self.exonIds = gtf_transcript.exonIds
        self.begin = gtf_transcript.begin
        self.end = gtf_transcript.end

    def add(self, gtf_transcript):
        """Adds Transcripts to GTF Gene"""

        self.gtf_transcripts.add(gtf_transcript)
        self.transcriptIds.add(gtf_transcript.id)
        self.exonIds |= gtf_transcript.exonIds
        self.begin = min(self.begin, gtf_transcript.begin)
        self.end = max(self.end, gtf_transcript.end)


def argument_parser():
    """Parses Agruments from Terminal"""

    global input_files
    global threshold

    parser = ArgumentParser()

    parser.add_argument("-r", dest="refgtf", help='''a set of known mRNAs to use as a reference for assessing 
       the accuracy of mRNAs or gene models given in <input.gtf>''',
                        metavar="<reference_mrna>", nargs='?', required=True)

    parser.add_argument("-i", help='''provide a text file with a list of GTF files to process instead
       of expecting them as command line arguments (useful when a large number
       of GTF files should be processed)''', metavar="<input_list>", nargs='?')

    parser.add_argument('input_annotations', help='''a set of known mRNAs to use as a reference for assessing 
       the accuracy of mRNAs or gene models given in <input.gtf>''', metavar='<input.gtf>', nargs='*')

    parser.add_argument('-th', help='a basepair threshold for exon matching',
                        metavar='basepair threshold', type=int, nargs=1, default=[0])

    args = parser.parse_args()
    if len(args.input_annotations) == 0 and args.i is None:
        print("Not enough arguments, enter at least one assembled annotation file")

    input_annotations = []

    for x in args.input_annotations:
        input_annotations.append(x)
    if args.i is not None:
        with open(args.i, 'r') as input_list:
            for line in input_list:
                input_annotations.append(line.rstrip())

    input_files = input_annotations
    reference = str(args.refgtf)
    threshold = args.th[0]
    return reference


def parse_database():
    """Converts GTF file to gffutills database files if necessary and parses all reference exons"""
    global input_files
    global threshold
    global reference_exons

    reference = argument_parser()

    # Converts reference GTF to DB if needed
    if not reference.endswith(".db"):
        gffutils.create_db(reference, dbfn=reference + '.db', force=True,
                           disable_infer_genes=True, disable_infer_transcripts=True)
        reference += '.db'

    # Converts query GTFs to DB if needed
    for x in range(len(input_files)):
        if not input_files[x].endswith(".db"):
            gffutils.create_db(input_files[x], dbfn=input_files[x] + '.db', force=True,
                               disable_infer_genes=True, disable_infer_transcripts=True)
            input_files[x] += '.db'

    # Parses reference DB
    refdb = gffutils.FeatureDB(reference)
    # Parses all reference exons
    reference_exons = refdb.features_of_type('exon')


def write_last_exon(cnt, dotexon, quer_intervals, last_exon):
    """Writes the last exon in .exon file"""
    # retrives the interval object to match it to reference
    for cinter in quer_intervals.search(last_exon.begin, last_exon.end, strict=True):
        if cinter.begin == last_exon.begin and cinter.end == last_exon.end:
            # get the best reference match
            bests = match_interval(cinter, last_exon.chrom, last_exon.strand)
            break
    # If a match (best match) was found write each match in .exon file
    if bests:
        for bintr, bval in bests.items():
            dotexon.write('{}\t{}\t{}-{}[{}]|{}[{}]\t{}\t{}-{}[{}]|{}\t{}\t({},{})\t({})\n'.format(
                cnt - 1, last_exon.chrom, cinter.begin, cinter.end - 1, cinter.data.strand, last_exon.transcript_id,
                cinter.data.transcriptIds[last_exon.transcript_id], bval[1], bintr.begin, bintr.end - 1,
                bintr.data.strand, '|'.join(['{}[{}]'.format(k, v) for k, v in bintr.data.transcriptIds.items()]),
                bval[0], bintr.begin - cinter.begin, cinter.end - bintr.end, NOTES[cinter.data.note]
            ))
    else:
        dotexon.write('{}\t{}\t{}-{}[{}]|{}[{}]\tNovel\t-\t-\t-\t-\n'.format(
            cnt - 1, last_exon.chrom, cinter.begin, cinter.end - 1, cinter.data.strand, last_exon.transcript_id,
            cinter.data.transcriptIds[last_exon.transcript_id]
        ))


def parse_reference():
    """Generates Reference stats and constructs reference dictionaries"""
    global reference_exons
    global ref_chrom_dict
    global ref_introns

    # Next line defines a variable which will be used to calculate exonic unique length
    ref_ulen = {'+': IntervalTree(), '-': IntervalTree(), '.': IntervalTree()}
    # Next line defines a variable which will be used to calculate reference statistics
    REF_STATS = {'Total_Exons': [], 'Unique_Exons': [], '5_Prime_UTR': 0, '3_Prime_UTR': 0, 'Intronic_Edges': 0}

    gene_strand_count = {}
    tran_strand_count = {}
    ref_chrom_dict = {}
    ref_gene_dict = {}
    ref_transcript_dict = {}
    number = 0
    intron_start = 0
    intron_end = 0
    last_interval = None
    prev_transcript_id = None

    # Loops over all reference exons
    for i, exon in enumerate(tqdm(reference_exons)):
        exon = Exon(exon)
        exon_interval = Interval(exon.begin, exon.end + 1)

        # Calculates Unique length
        mn = exon.begin
        mx = exon.end + 1
        # Loops over all overlapping intervals removes them
        for b, e, d in ref_ulen[exon.strand][mn - 1:mx]:
            if b < mn:
                mn = b
            if e > mx:
                mx = e
            ref_ulen[exon.strand].removei(b, e, d)
        # Adds and interval to unique length starting from smallest start to the largest end
        ref_ulen[exon.strand].addi(mn, mx, i)

        # Checks if the current exon in a new transcript
        if prev_transcript_id != exon.transcript_id:
            if last_interval:
                if number == 'First':
                    last_interval.transcriptIds[prev_transcript_id] = 'Single'
                else:
                    last_interval.transcriptIds[prev_transcript_id] = 'Last'

                REF_STATS['Intronic_Edges'] -= 1
                REF_STATS['3_Prime_UTR'] += 1
            number = 'First'
            prev_transcript_id = exon.transcript_id

            if exon.strand == '+':
                intron_start = exon.end
                intron_end = None
            else:
                intron_start = None
                intron_end = exon.start
        else:
            if number == 'First':
                number = 1
            number += 1

            if exon.strand == '-':
                intron_start = exon.end
            else:
                intron_end = exon.start

        REF_STATS['Total_Exons'].append(exon.size)
        if number == 'Single':
            REF_STATS['5_Prime_UTR'] += 1
            REF_STATS['3_Prime_UTR'] += 1
        elif number == 'First':
            REF_STATS['5_Prime_UTR'] += 1
            REF_STATS['Intronic_Edges'] += 1
        elif number == 'Last':
            REF_STATS['Intronic_Edges'] += 1
            REF_STATS['3_Prime_UTR'] += 1
        else:
            REF_STATS['Intronic_Edges'] += 1
            REF_STATS['Intronic_Edges'] += 1

        if intron_start and intron_end:
            ref_introns[intron_start:intron_end] = exon.strand
            if exon.strand == '-':
                intron_end = exon.start
            else:
                intron_start = exon.end

        # To Replace gene id with gene name
        if exon.attributes.get('gene_name'):
            gene_name = exon.attributes.get('gene_name')[0]
            gene_name_id[exon.gene_id] = gene_name

        # Creates a dictionary for the chromosome if it didn't exist
        if exon.chrom not in ref_chrom_dict:
            ref_chrom_dict[exon.chrom] = {}

        # Creates a dictionary for the strand if it didn't exist
        if exon.strand not in ref_chrom_dict[exon.chrom]:
            # Creates an interval tree for intervals and exons dictionary relative to the current chromosome and strand
            ref_intervals = IntervalTree()
            ref_exons_dict = {}

            last_interval = GtfInterval(exon_interval, exon, number)
            ref_intervals[exon.start: exon.end + 1] = last_interval
            # REF_STATS['Unique_Exons'].append(exon.size)

            last_exon = GtfExon(exon.id, last_interval, exon.transcript_id, exon.gene_id)
            ref_exons_dict[exon.id] = last_exon

        else:
            # Retrives the interval tree and exon dictionary
            ref_intervals, ref_exons_dict = ref_chrom_dict[exon.chrom][exon.strand]
            found = False
            dz = None
            # Checks if the interval already exists
            for b, e, z in ref_intervals.search(exon.begin, exon.end + 1, strict=True):
                if b == exon_interval.begin and e == exon_interval.end:
                    found = True
                    dz = z
                    break
            # Means the interval didn't exist
            if not found:
                last_interval = GtfInterval(exon_interval, exon, number)
                ref_intervals[exon.begin: exon.end + 1] = last_interval
                # REF_STATS['Unique_Exons'].append(exon.size)
            else:
                if dz is None:
                    print('Interval not found')
                    exit(1)
                # Add exon data to the existing interval
                dz.add(exon.id, exon.transcript_id, number, exon.gene_id)
                last_interval = dz

            if exon.id not in ref_exons_dict:
                last_exon = GtfExon(exon.id, last_interval, exon.transcript_id, exon.gene_id)
                ref_exons_dict[exon.id] = last_exon
            else:
                last_exon = ref_exons_dict[exon.id]
                print('Duplicate Exon ID', exon.id, exon_interval)

        ref_chrom_dict[exon.chrom][exon.strand] = [ref_intervals, ref_exons_dict]

        if exon.transcript_id not in ref_transcript_dict:
            last_transcript = GtfTranscript(exon.transcript_id, last_exon, exon.gene_id)
            ref_transcript_dict[exon.transcript_id] = last_transcript
        else:
            last_transcript = ref_transcript_dict[exon.transcript_id]
            if last_transcript.gene_id != exon.gene_id:
                print('Duplicate Transcript id in different genes')
            elif last_transcript.id != exon.transcript_id:
                print('Duplicate Transcript', last_transcript.id)
            last_transcript.add(last_exon)

        if exon.gene_id not in ref_gene_dict:
            ref_gene_dict[exon.gene_id] = GtfGene(exon.gene_id, last_transcript)
        else:
            if ref_gene_dict[exon.gene_id].id != exon.gene_id:
                print('Duplicate Gene')
            ref_gene_dict[exon.gene_id].add(last_transcript)

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
    else:
        last_interval.transcriptIds[prev_transcript_id] = 'Last'

    REF_STATS['Intronic_Edges'] -= 1
    REF_STATS['3_Prime_UTR'] += 1

    for chromo in ref_chrom_dict:
        for strand in ref_chrom_dict[chromo]:
            ref_chrom_dict[chromo][strand].append({})
            ref_chrom_dict[chromo][strand].append({})
    # Determines each transcript's strand
    for chromo in tran_strand_count:
        for tran, strands in tran_strand_count[chromo].items():
            if strands['+'] > strands['.'] and strands['+'] > strands['-']:
                strand = '+'
            elif strands['-'] > strands['.'] and strands['-'] > strands['+']:
                strand = '-'
            else:
                strand = '.'
            ref_chrom_dict[chromo][strand][2][tran] = ref_transcript_dict[tran]
    # Determines each gene's strand
    for chromo in gene_strand_count:
        for gene, strands in gene_strand_count[chromo].items():
            if strands['+'] > strands['.'] and strands['+'] > strands['-']:
                strand = '+'
            elif strands['-'] > strands['.'] and strands['-'] > strands['+']:
                strand = '-'
            else:
                strand = '.'
            ref_chrom_dict[chromo][strand][3][gene] = ref_gene_dict[gene]

    # Write temporarily Reference stats file
    with open('refstat.txt', 'w') as refstatf:
        wf = False
        for st in ['+', '-', '.']:
            for b, e, d in ref_ulen[st]:
                REF_STATS['Unique_Exons'].append(e - b)
        for key, value in REF_STATS.items():
            if 'Exons' in key:
                val = np.array(value)
                refstatf.write('{}(count|total|mean|std)\t{}|{}|{}|{}\n'.
                               format(key, len(val), val.sum(), round(val.mean(), 1), round(val.std(), 1)))
            else:
                if key == '5_Prime_UTR':
                    if value == REF_STATS['3_Prime_UTR']:
                        continue
                    else:
                        print('Something went wrong with stats UTRS don\'t match')
                        wf = True
                if key == '3_Prime_UTR' and not wf:
                    key = 'Transcripts'
                refstatf.write('{}(count)\t{}\n'.format(key, value))
    # return ref_chrom_dict


def parse_query(num):
    """Constructs query dictionaries and generates .exon files
       num is the number of the current query file
    """
    global interval_matches
    global interval_best_matches
    global ref_matches
    global ref_best_matches

    # Resets matches dictionaries
    interval_matches = {}
    interval_best_matches = {}
    ref_matches = {}
    ref_best_matches = {}

    # Parses all query exons
    query_exons = gffutils.FeatureDB(input_files[num]).features_of_type('exon')

    # Creates .exon file and writes the header
    dotexon = open(input_files[num].split('/')[-1] + '.exon', 'w')
    dotexon.write("ExonID\tChromosome\tQuery(Coordinates[strand]|Transcript[exon_number])\tMatch_Type\t" +
                  "Reference(Best_Match_Coordinates|Transcript[exon_number])\tShared\tBase_Difference\tNotes\n")

    gene_strand_count = {}
    tran_strand_count = {}
    quer_chrom_dict = {}
    quer_transcript_dict = {}
    quer_gene_dict = {}
    number = 0
    last_interval = None
    prev_transcript_id = None
    prev_strand = None

    # Loops over all query exons
    for cnt, exon in enumerate(tqdm(query_exons)):
        cnt += 1
        exon = Exon(exon)
        exon_interval = Interval(exon.begin, exon.end + 1)

        # Checks if the current exon in a new transcript
        if prev_transcript_id != exon.transcript_id:
            if last_interval:
                if number == 'First' or (number == 'Last' and prev_strand == '-'):
                    last_interval.transcriptIds[prev_transcript_id] = 'Single'
                else:
                    if prev_strand == '-' and reverse:
                        last_interval.transcriptIds[prev_transcript_id] = 'First'
                    else:
                        last_interval.transcriptIds[prev_transcript_id] = 'Last'
            if exon.strand == '-':
                number = 'Last'
                assumed = True
                reverse = True
            else:
                number = 'First'
            prev_transcript_id = exon.transcript_id
        else:
            if exon.strand == '-':
                if assumed:
                    if prev_strand == '-' and last_interval.end > exon.begin:
                        last_interval.transcriptIds[exon.transcript_id] = 'First'
                        reverse = False
                    assumed = False
            number = 'Mid'

        if cnt > 1:
            write_last_exon(cnt, dotexon, quer_intervals, last_exon)

        # Creates a dictionary for the chromosome if it didn't exist
        if exon.chrom not in quer_chrom_dict:
            quer_chrom_dict[exon.chrom] = {}

        # Creates a dictionary for the strand if it didn't exist
        if exon.strand not in quer_chrom_dict[exon.chrom]:
            # Creates an interval tree for intervals and exons dictionary relative to the current chromosome and strand
            quer_intervals = IntervalTree()
            quer_exons_dict = {}

            last_interval = GtfInterval(exon_interval, exon, number)
            quer_intervals[exon.start: exon.end + 1] = last_interval

            last_exon = GtfExon(exon.id, last_interval, exon.transcript_id, exon.gene_id)
            quer_exons_dict[exon.id] = last_exon

        else:
            # Retrives the interval tree and exon dictionary
            quer_intervals, quer_exons_dict = quer_chrom_dict[exon.chrom][exon.strand]
            found = False
            dz = None
            # Checks if the interval already exists
            for b, e, z in quer_intervals.search(exon.begin, exon.end + 1, strict=True):
                if b == exon_interval.begin and e == exon_interval.end:
                    found = True
                    dz = z
                    break
            # Means the interval didn't exist
            if not found:
                last_interval = GtfInterval(exon_interval, exon, number)
                quer_intervals[exon.begin: exon.end + 1] = last_interval
            else:
                if dz is None:
                    print('Interval not found')
                    exit(1)
                # Add exon data to the existing interval
                dz.add(exon.id, exon.transcript_id, number, exon.gene_id)
                last_interval = dz

            if exon.id not in quer_exons_dict:
                last_exon = GtfExon(exon.id, last_interval, exon.transcript_id, exon.gene_id)
                quer_exons_dict[exon.id] = last_exon
            else:
                last_exon = quer_exons_dict[exon.id]
                print('Duplicate Exon ID', exon.id, exon_interval)

        quer_chrom_dict[exon.chrom][exon.strand] = [quer_intervals, quer_exons_dict]

        if exon.transcript_id not in quer_transcript_dict:
            last_transcript = GtfTranscript(exon.transcript_id, last_exon, exon.gene_id)
            quer_transcript_dict[exon.transcript_id] = last_transcript
        else:
            last_transcript = quer_transcript_dict[exon.transcript_id]
            if last_transcript.gene_id != exon.gene_id:
                print('Duplicate Transcript id in different genes')
            elif last_transcript.id != exon.transcript_id:
                print('Duplicate Transcript', last_transcript.id)
            last_transcript.add(last_exon)

        if exon.gene_id not in quer_gene_dict:
            quer_gene_dict[exon.gene_id] = GtfGene(exon.gene_id, last_transcript)
        else:
            if quer_gene_dict[exon.gene_id].id != exon.gene_id:
                print('Duplicate Gene')
            quer_gene_dict[exon.gene_id].add(last_transcript)

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

        prev_strand = exon.strand

    if number == 'First' or (number == 'Last' and prev_strand == '-'):
        last_interval.transcriptIds[prev_transcript_id] = 'Single'
    else:
        if prev_strand == '-' and reverse:
            last_interval.transcriptIds[prev_transcript_id] = 'First'
        else:
            last_interval.transcriptIds[prev_transcript_id] = 'Last'

    for chromo in quer_chrom_dict:
        for strand in quer_chrom_dict[chromo]:
            quer_chrom_dict[chromo][strand].append({})
            quer_chrom_dict[chromo][strand].append({})
    # Determines each transcript's strand
    for chromo in tran_strand_count:
        for tran, strands in tran_strand_count[chromo].items():
            if strands['+'] > strands['.'] and strands['+'] > strands['-']:
                strand = '+'
            elif strands['-'] > strands['.'] and strands['-'] > strands['+']:
                strand = '-'
            else:
                strand = '.'
            quer_chrom_dict[chromo][strand][2][tran] = quer_transcript_dict[tran]
    # Determines each gene's strand
    for chromo in gene_strand_count:
        for gene, strands in gene_strand_count[chromo].items():
            if strands['+'] > strands['.'] and strands['+'] > strands['-']:
                strand = '+'
            elif strands['-'] > strands['.'] and strands['-'] > strands['+']:
                strand = '-'
            else:
                strand = '.'
            quer_chrom_dict[chromo][strand][3][gene] = quer_gene_dict[gene]
    # Write the very last exon to .exon file
    write_last_exon(cnt + 1, dotexon, quer_intervals, last_exon)
    return quer_chrom_dict


def is_overlapped(feature_1, feature_2):
    """Calculates the size of overlap and judgement of that overlap if exists"""

    judgement = 'No_Overlap'
    if type(feature_1) is GtfExon:
        start_1 = feature_1.interval.begin
        end_1 = feature_1.interval.end
    elif type(feature_1) is Interval:
        start_1 = feature_1.begin
        end_1 = feature_1.end
    else:
        start_1 = feature_1.start
        end_1 = feature_1.end
    if type(feature_2) is GtfExon:
        start_2 = feature_2.interval.begin
        end_2 = feature_2.interval.end
    elif type(feature_2) is Interval:
        start_2 = feature_2.begin
        end_2 = feature_2.end
    else:
        start_2 = feature_2.start
        end_2 = feature_2.end
    size_1 = end_1 - start_1
    size_2 = end_2 - start_2
    if start_1 > start_2 or (start_1 == start_2 and size_2 > size_1):
        temp_s, temp_e, temp_z = start_1, end_1, size_1
        start_1, end_1, size_1 = start_2, end_2, size_2
        start_2, end_2, size_2 = temp_s, temp_e, temp_z
    if start_2 > end_1:
        return [None, judgement]
    if start_1 == start_2 and end_1 == end_2:
        judgement = 'Pefect_Match'
    else:
        if end_2 > end_1:
            judgement = 'Overlapped'
        else:
            judgement = "Contained"
    if end_2 > end_1:
        return [end_1 - start_2, judgement]
    return [size_2, judgement]


# def match_isoforms(ref_isos, query_iso, query_chrom_dict, chrom, strand):
#     global ref_matches
#     ref_iso_ids = dict()
#     ref_iso_intervals = dict()
#     for ref_iso in ref_isos:
#         ref_iso_ids[ref_iso] = ref_chrom_dict[chrom][strand][2][ref_iso].exonIds
#         ref_iso_intervals[ref_iso] = set()
#         for id in ref_iso_ids:
#             ref_iso_intervals[ref_iso] |= set(ref_chrom_dict[chrom][strand][1][id].interval)
#     query_iso_ids = sorted(query_iso.exonIds)
#     query_iso_intervals = [query_chrom_dict[chrom][strand][1][id].interval for id in query_iso_ids]
#     query_iso_intervals = set(query_iso_intervals)
#     overlaps = [ref_matches.get(intr, None) for intr in query_iso_intervals]
#     print(overlaps, '\nafter\n')
#     overlaps = [v for x in overlaps for k, v in x.items() if v in ref_iso_intervals]
#     print(overlaps)
#     exit()


# def find_super_loci(quer_chrom_dict):
#     global ref_chrom_dict
#     loc_ind = {}
#     ind = 0
#     linked = {}
#     loci = {}
#     super_loci = {}
#     quer_loci_gene = {}
#     for chrom in quer_chrom_dict:
#         for strand in quer_chrom_dict[chrom]:
#             quer_intervals, quer_exons_dict, quer_transcript_dict, quer_gene_dict =\
#                 quer_chrom_dict[chrom][strand]
#             for locus in quer_gene_dict:
#                 quer = locus
#                 locus = quer_gene_dict[locus]
#                 interval = Interval(locus.start, locus.end - 1)
#                 found_overlap = False
#                 if chrom in ref_chrom_dict:
#                     if strand in ref_chrom_dict[chrom]:
#                     # for st in ref_chrom_dict[chrom]:
#                         ref = get_all(ref_chrom_dict[chrom][strand], 'genes', interval)
#                         if ref:
#                             # if st != '.' and strand != '.':
#                             #     if st != strand:
#                             #         print('found a locus overlap on opposite strand reference strand is ' + st \
#                             #               + ' interval is ' + str(interval) + ' ' + strand + '\n', ref)
#                             found_overlap = True
#                             for loc in ref:
#                                 if loc in loc_ind:
#                                     prev_ind = loc_ind[loc]
#                                     break
#                                 else:
#                                     prev_ind = 0
#                             if prev_ind == 0:
#                                 ind += 1
#                                 linked[ind] = [set(), set()]
#                             for loc in ref:
#                                 if prev_ind != 0:
#                                     loc_ind[loc] = prev_ind
#                                     linked[prev_ind][1].add(loc)
#                                     linked[prev_ind][0].add(quer)
#                                 else:
#                                     loc_ind[loc] = ind
#                                     linked[ind][1].add(loc)
#                                     linked[ind][0].add(quer)
#
#                             if quer not in quer_loci_gene:
#                                 quer_loci_gene[quer] = [chrom, strand]
#
#                 if not found_overlap:
#                     if quer not in loci:
#                         loci[quer] = [chrom, strand, interval]
#                     else:
#                         if loci[quer][0] != chrom or loci[quer][1] != strand or loci[quer][2] != interval:
#                             print("Something went wrong dual gene", quer, interval)
#                             exit(1)
#     for v in linked.values():
#         v[0] = list(v[0])
#         if v[0]:
#             # print(v, ' v\nv0 ', v[0], '\nv1 ', v[1])
#             super_loci[tuple(v[0])] = quer_loci_gene[v[0][0]] + list(v[1])
#     return [loci, super_loci]


# def write_loci_file(quer_chrom_dict, loci, super_loci, x):
#     formatted = {}
#     super_loci_counter = 1
#     with open('gtfcomp.' + input_files[x].split('/')[-1] + '.loci', 'w') as f:
#         for k, v in sorted(loci.items()):
#             chrom, strand, interval = v
#             line = k+'\t'+chrom+'['+strand+']'+str(interval.begin)+'-'+str(interval.end)+'\t-\t'
#             # f.write(line)
#             quer_gene_dict = quer_chrom_dict[chrom][strand][3]
#             quer = sorted(quer_gene_dict[k].transcriptIds, reverse=True)
#             for l in range(len(quer) - 1):
#                 line += quer.pop() + ','
#                 # f.write(quer.pop() + ',')
#             # f.write(quer.pop() + '\n')
#             line += quer.pop() + '\n'
#             if chrom not in formatted:
#                 formatted[chrom] = {}
#             if strand not in formatted[chrom]:
#                 formatted[chrom][strand] = {}
#             if interval.begin not in formatted[chrom][strand]:
#                 formatted[chrom][strand][interval.begin] = [line]
#             else:
#                 formatted[chrom][strand][interval.begin].append(line)
#         for keys, v in super_loci.items():
#             quer = set()
#             chrom, strand, genes = v[0], v[1], v[2:]
#             quer_gene_dict = quer_chrom_dict[chrom][strand][3]
#             ref = ''
#             start = 1000**1000
#             end = 0
#             for k in keys:
#                 if k in quer_gene_dict:
#                     gen = quer_gene_dict[k]
#                     quer |= gen.transcriptIds
#                     start = min(start, gen.start)
#                     end = max(end, gen.end)
#             if quer:
#                 quer = sorted(quer, reverse=True)
#                 k = 'SLOC_'+str(super_loci_counter).zfill(6)
#                 super_loci_counter += 1
#                 ref_gene_dict = ref_chrom_dict[chrom][strand][3]
#                 for g in genes:
#                     if g in ref_gene_dict:
#                         gen = ref_gene_dict[g]
#                         g = gene_name_id.get(g, g)
#                         start = min(start, gen.start)
#                         end = max(end, gen.end)
#                         ref += g + '|'
#                         trans = sorted(gen.transcriptIds, reverse=True)
#                         for l in range(len(trans) - 1):
#                             ref += trans.pop() + '|'
#                         ref += trans.pop() + ','
#                 if ref:
#                     ref = ref[:-1]
#                     line = k+'\t'+chrom+'['+strand+']'+str(start)+'-'+str(end - 1)+'\t'+ref+'\t'
#                     # f.write(line)
#                     for l in range(len(quer) - 1):
#                         line += quer.pop() + ','
#                         # f.write(quer.pop() + ',')
#                     # f.write(quer.pop() + '\n')
#                     line += quer.pop() + '\n'
#                     if chrom not in formatted:
#                         formatted[chrom] = {}
#                     if strand not in formatted[chrom]:
#                         formatted[chrom][strand] = {}
#                     if start not in formatted[chrom][strand]:
#                         formatted[chrom][strand][start] = [line]
#                     else:
#                         formatted[chrom][strand][start].append(line)
#         for ch in sorted(formatted):
#             for st in sorted(formatted[ch]):
#                 for intr in sorted(formatted[ch][st]):
#                     for item in formatted[ch][st][intr]:
#                         f.write(item)


def match_interval(interval, chrom, stran):
    """Matches query interval to the reference"""

    global interval_matches
    global interval_best_matches
    global ref_matches
    global ref_best_matches
    global ref_chrom_dict

    # If interval was previously matched return the match
    if interval in interval_matches:
        return interval_best_matches[interval]

    # Indicates if any match was found
    found = False

    # No note initially
    note = 0

    if chrom in ref_chrom_dict:
        if stran in ref_chrom_dict[chrom]:
            # Search Same strand first
            found = ref_chrom_dict[chrom][stran][0][interval]

        # If nothing was found search unknown first then opposite strand (strand is either + or -)
        if not found and stran != '.':
            # Search unknown first
            if '.' in ref_chrom_dict[chrom]:
                found = ref_chrom_dict[chrom]['.'][0][interval]
                # Note 1: A match was found on the unknown strand of the reference
                note = 1
            # Then search opposite
            if not found:
                if stran == '+' and '-' in ref_chrom_dict[chrom]:
                    found = ref_chrom_dict[chrom]['-'][0][interval]

                elif stran == '-' and '+' in ref_chrom_dict[chrom]:
                    found = ref_chrom_dict[chrom]['+'][0][interval]

                if found:
                    # Note 2: A match was found on the opposite strand on the reference
                    note = 2
                    # print('found match for ', interval, ' on reverse strand on reference')

        # If strand was unknown search both strands
        elif not found and stran == '.':
            found_pos = set()
            found_neg = set()

            if '+' in ref_chrom_dict[chrom]:
                found_pos = ref_chrom_dict[chrom]['+'][0][interval]

            if '-' in ref_chrom_dict[chrom]:
                found_neg = ref_chrom_dict[chrom]['-'][0][interval]

            if found_pos and found_neg:
                # print('unknown of both strands')
                # Note 3: A match was found on the + and - strands of the reference query's strand is unknown
                note = 3

            found = found_pos | found_neg

        if found:
            # Convert the found set to a list
            comp = [x for x in found]

            # Compute overlaps and judgments of each match (reference match)
            comp_overlaps = list(map(is_overlapped, comp, itertools.repeat(interval, len(comp))))

            # Unpacks the previous step 2D list into 2 1D lists
            comp_overlaps, comp_judgements = zip(*comp_overlaps)
            # comp_judgements = [x[1] for x in comp_overlaps]
            # comp_overlaps = [x[0] for x in comp_overlaps]

            # Computes the size of each match (reference match)
            comp_sizes = list(map(lambda xx: xx.end - xx.begin, comp))

            # Picks the best match from all matches
            ii, ovlp = [0], 0
            for i, v in enumerate(comp_overlaps):
                if v > ovlp or (v == ovlp and comp_sizes[i] < comp_sizes[ii[0]]):
                    ovlp = v
                    ii = [i]
                elif v == ovlp and comp_sizes[i] == comp_sizes[ii[0]]:
                    ii.append(i)

            # Query interval size
            intr_size = interval.end - 1 - interval.begin
            if len(ii) > 1:
                if note == 0:
                    # Note 4: More than one best match was found all on the same strand same as query strand
                    note = 4
                elif note == 1:
                    # Note 5: More than one best match was found on the unknown strand of reference
                    note = 5
                elif note == 2:
                    # Note 6: More than one best match was found all on the same strand opposite of query strand
                    note = 6
                elif note == 3:
                    # Note 7: More than one best match was found on both strands of reference query's strand is unknown
                    note = 7

            # Sets the note for the interval if found
            interval.data.note = note

            # Checks if something went wrong and matched the interval twice
            # Saves the best match in query best matches dictionary
            if interval not in interval_best_matches:
                interval_best_matches[interval] = {}
                for xi in ii:
                    interval_best_matches[interval][comp[xi]] = (ovlp, comp_judgements[xi])
            else:
                print('Dual Interval found ', interval, 'more than once ', interval_best_matches[interval],
                      '\n', (comp[ii[0]], ovlp, comp_sizes[ii[0]], comp_judgements[ii[0]]))

            # To Calculate best match for the refrenece and to save all possible matches
            for i, v in enumerate(comp):
                if v in ref_best_matches:
                    ri = ref_best_matches[v]
                    rii = six.next(six.itervalues(ri))
                    if (comp_overlaps[i] > rii[0]) or (comp_overlaps[i] == rii[0] and intr_size < rii[2]):
                        ri = {interval: (comp_overlaps[i], comp_judgements[i], intr_size)}
                        # ri[interval] = (comp_overlaps[i], comp_judgements[i])
                    elif comp_overlaps[i] == rii[0] and comp_sizes[i] == rii[1]:
                        ri[interval] = (comp_overlaps[i], comp_judgements[i], intr_size)
                else:
                    ref_best_matches[v] = {}
                    ri = ref_best_matches[v]
                    ri[interval] = (comp_overlaps[i], comp_judgements[i], intr_size)

                if v not in ref_matches:
                    ref_matches[v] = {}
                ref_matches[v][interval] = (comp_overlaps[i], comp_judgements[i])

                if interval not in interval_matches:
                    interval_matches[interval] = {}
                interval_matches[interval][v] = (comp_overlaps[i], comp_judgements[i])

            return interval_best_matches[interval]


def write_itrack(num, quer_chrom_dict):
    """Writes .itracking files
       num is the number of the current query file"""

    global interval_matches
    global interval_best_matches
    global ref_matches
    global ref_best_matches
    global ref_introns

    # Create .itracking file
    itrack = open(input_files[num].split('/')[-1] + '.itracking', 'w')
    # Write header
    itrack.write("ID\tChromosome\tQuery(Coordinates[strand]|Transcript[exon_number])\t"
                 "Reference(Ref_Coordinates[strand]|Transcript[exon_number])\tType\t(Start,End)\t[Query_Structure]"
                 "\t[Reference_Structure]\tNotes\n")
    # Index for each row
    cnt = 0
    # To hold reported intervals not to report them again
    reported = set()
    for chrom in quer_chrom_dict:
        for strand in ['.', '+', '-']:
            if strand not in quer_chrom_dict[chrom]:
                continue
            for interval in sorted(quer_chrom_dict[chrom][strand][0]):
                if interval in reported:
                    continue
                cnt += 1
                if interval in interval_matches:
                    qset = {interval}
                    reported.add(interval)
                    rset = set()
                    # Get the first reference best match's strand
                    ref_st = six.next(six.iterkeys(interval_best_matches[interval])).data.strand

                    if interval.data.strand != strand:
                        print('SOMETHING IS WRONG STRANDS DOESN\'t MATCH Code: 1')

                    # Get current strand
                    quer_st = interval.data.strand

                    # Add note for the unknown strand if found on both reference strands
                    if len(interval_best_matches[interval]) > 1 and interval.data.note == 7:
                        note = '(A best match on both strands was found only used {} (arbitrary decision))'.format(
                            ref_st)
                    else:
                        note = '({})'.format(NOTES[interval.data.note])

                    # Loop over the query intervals and find all reference matches and query matches
                    tempset = qset.copy()
                    tempset2 = True
                    while tempset2:
                        if not tempset:
                            tempset = tempset2
                        tempset2 = set()
                        for qinterval in tempset:
                            for rinterval in interval_matches[qinterval]:
                                # Skip if reference strand doesn't match best match strand
                                if rinterval.data.strand != ref_st and rinterval.data.strand != '.':
                                    continue
                                if rinterval not in rset:
                                    rset.add(rinterval)
                                    for qrinterval in ref_matches[rinterval]:
                                        # Skip if already reported
                                        if qrinterval in reported:
                                            continue
                                        # Skip if query strand doesn't match the current strand
                                        if qrinterval.data.strand != quer_st and qinterval.data.strand != '.':
                                            continue
                                        # Skip if a match was found on opposite strand
                                        if qrinterval.data.note == 2:
                                            continue
                                        if qrinterval not in qset:
                                            tempset2.add(qrinterval)
                                            qset.add(qrinterval)
                                            # Add current interval to reported set
                                            reported.add(qrinterval)
                        tempset = set()
                    qset = sorted(qset)
                    rset = sorted(rset)
                    cord = ''
                    qstat = '-'
                    rstat = '-'
                    if len(qset) == 1:
                        mtype = 'SQ'
                        qm = next(iter(qset))
                        cord = '{}-{}[{}]|{}'.format(
                                qm.begin, qm.end - 1, qm.data.strand,
                                '|'.join(['{}[{}]'.format(k, v) for k, v in qm.data.transcriptIds.items()]))
                        qstart = qm.begin
                        qend = qm.end

                    elif len(qset) > 1:
                        mtype = 'MQ'
                        qstat = '['
                        for i, qm in enumerate(qset):
                            cord += '{}-{}[{}]|{}, '.format(
                                qm.begin, qm.end - 1, qm.data.strand,
                                '|'.join(['{}[{}]'.format(k, v) for k, v in qm.data.transcriptIds.items()]))
                            if i > 0:
                                qstat += str(prev_end - qm.begin + 1) + " "
                            else:
                                qstart = qm.begin
                                qend = qm.end
                            if qm.end <= qend and i != 0:
                                if note == '(-)':
                                    note = '(Total Overlap occured at query number '
                                elif note[-1] == ')':
                                    note = note[:-1] + ' Total Overlap occured at query number '
                                note += '{}, '.format(i)
                            qstart = min(qstart, qm.begin)
                            qend = max(qend, qm.end)
                            prev_end = qm.end - 1

                        cord = cord[:-2]
                        qstat = qstat[:-1] + ']'
                        note = note[:-2] + ')' if note[-1] == ' ' else note
                    cord += '\t'

                    if len(rset) == 1:
                        mtype += 'SR'
                        rm = next(iter(rset))
                        cord += 'r|{}-{}[{}]|{}'.format(
                            rm.begin, rm.end - 1, rm.data.strand,
                            '|'.join(['{}[{}]'.format(k, v) for k, v in rm.data.transcriptIds.items()]))
                        rstart = rm.begin
                        rend = rm.end

                    elif len(rset) > 1:
                        mtype += 'MR'
                        rstat = '['
                        for i, rm in enumerate(rset):
                            cord += 'r|{}-{}[{}]|{}, '.format(
                                rm.begin, rm.end - 1, rm.data.strand,
                                '|'.join(['{}[{}]'.format(k, v) for k, v in rm.data.transcriptIds.items()]))
                            if i > 0:
                                rstat += str(prev_end - rm.begin + 1) + " "
                            else:
                                rstart = rm.begin
                                rend = rm.end
                            if rm.end <= rend and i != 0:
                                if note == '(-)':
                                    note = '(Total Overlap occured at reference number '
                                elif note[-1] == ')':
                                    note = note[:-1] + ' Total Overlap occured at reference number '
                                note += '{}, '.format(i)
                            rstart = min(rstart, qm.begin)
                            rend = max(rend, rm.end)
                            prev_end = rm.end - 1

                        cord = cord[:-2]
                        rstat = rstat[:-1] + ']'
                        note = note[:-2] + ')' if note[-1] == ' ' else note

                    if mtype == 'SQSR':
                        if note != '-':
                            note = '({}, {})'.format(note.replace('(', '').replace(')', ''), NOTES[qm.data.note])
                        elif qm.data.note != 0:
                            note = '({})'.format(NOTES[qm.data.note])

                    start = rstart - qstart
                    end = qend - rend
                    itrack.write("id_{}\t{}\t{}\t".format(cnt, chrom, cord))
                    itrack.write("{}\t({},{})\t{}\t{}\t{}\n".format(mtype, start, end, qstat, rstat, note))
                else:
                    lvl = 'Intergenic'
                    for b, e, s in ref_introns[interval]:
                        if s == strand or strand == '.':
                            lvl = 'Intronic'
                    if interval.data.strand != strand:
                        print('SOMETHING IS WRONG STRANDS DOESN\'t MATCH Code: 2')
                    itrack.write("id_{}\t{}\t{}-{}[{}]|{}\t-\t".format(
                        cnt, chrom, interval.begin, interval.end - 1, interval.data.strand,
                        '|'.join(['{}[{}]'.format(k, v) for k, v in interval.data.transcriptIds.items()])
                    ))
                    itrack.write("SQNR\t-\t-\t-\t[{}]\n".format(lvl))


argument_parser()
parse_database()
parse_reference()

for x in range(len(input_files)):
    print(input_files[x])
    query_dict = parse_query(x)
    write_itrack(x, query_dict)
    # loci, super_loci = find_super_loci(quer_chrom_dict)
    # write_loci_file(quer_chrom_dict, loci, super_loci, x)
print("total time ", time.time() - time_begin)
