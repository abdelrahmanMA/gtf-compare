import gffutils
from argparse import ArgumentParser
from intervaltree import Interval, IntervalTree
import time
from tqdm import tqdm
import itertools
import operator

input_files = None


class Exon(gffutils.Feature):

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


class GtfInterval:

    def __init__(self, exon_id, transcript_id, gene_id, number):

        self.exonIds = set([exon_id])
        self.transcriptIds = set([transcript_id])
        self.geneIds = set([gene_id])
        self.tran_number = str(transcript_id) + "[" + str(number) + "]"

    def add(self, exon_id, transcript_id, gene_id, number):

        if exon_id not in self.exonIds:
            self.exonIds.add(exon_id)

        if transcript_id not in self.transcriptIds:
            self.transcriptIds.add(transcript_id)
            self.tran_number += "|" + str(transcript_id) + "[" + str(number) + "]"

        if gene_id not in self.geneIds:
            self.geneIds.add(gene_id)


class GtfExon:

    def __init__(self, interval, transcript_id, gene_id, number):

        self.interval = interval
        self.transcriptIds = set([transcript_id])
        self.geneIds = set([gene_id])
        self.begin = interval.begin
        self.end = interval.end

    def add(self, transcript_id, gene_id, number):

        if transcript_id not in self.transcriptIds:
            self.transcriptIds.add(transcript_id)

        if gene_id not in self.geneIds:
            self.geneIds.add(gene_id)


class GtfTranscript:

    def __init__(self, interval, exon_id, gene_id):

        self.intervals = IntervalTree()
        self.exonIds = set([exon_id])
        self.geneIds = set([gene_id])
        self.intervals.add(interval)
        self.start = interval.begin
        self.end = interval.end

    def add(self, interval, exon_id, gene_id):

        self.intervals.add(interval)
        self.start = min(self.start, interval.begin)
        self.end = max(self.end, interval.end)

        if exon_id not in self.exonIds:
            self.exonIds.add(exon_id)

        if gene_id not in self.geneIds:
            self.geneIds.add(gene_id)


class GtfGene:

    def __init__(self, interval, exon_id, transcript_id):

        self.intervals = IntervalTree()
        self.exonIds = set([exon_id])
        self.transcriptIds = set([transcript_id])
        self.intervals.add(interval)
        self.start = interval.begin
        self.end = interval.end

    def add(self, interval, exon_id, transcript_id):

        self.intervals.add(interval)
        self.start = min(self.start, interval.begin)
        self.end = max(self.end, interval.end)

        if exon_id not in self.exonIds:
            self.exonIds.add(exon_id)

        if transcript_id not in self.transcriptIds:
            self.transcriptIds.add(transcript_id)


def argument_parser():
    global time1
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

    parser.add_argument('-it', help='useful to speed up database generation if your GTF file has transcript features',
                        action="store_true", default=False)
    parser.add_argument('-ig', help='useful to speed up database generation if your GTF file has gene features',
                        action="store_true", default=False)

    args = parser.parse_args()
    if len(args.input_annotations) == 0 and args.i is None:
        print "Not enough arguments, enter at least one assembled annotation file"

    input_annotations = []

    for x in args.input_annotations:
        input_annotations.append(x)
    if args.i is not None:
        with open(args.i, 'r') as input_list:
            for line in input_list:
                input_annotations.append(line.rstrip())

    time1 = time.time()
    input_files = input_annotations
    reference = str(args.refgtf)
    threshold = args.th[0]
    op_infer_transcripts = args.it
    op_infer_genes = args.ig
    return [reference, op_infer_transcripts, op_infer_genes]


def parse_database():
    global input_files
    global threshold
    global reference_exons
    global query_exons
    reference, op_infer_transcripts, op_infer_genes = argument_parser()

    if not reference.endswith(".db"):
        gffutils.create_db(reference, dbfn=reference + '.db', force=True,
                           disable_infer_genes=op_infer_genes, disable_infer_transcripts=op_infer_transcripts)
        reference += '.db'

    for x in range(len(input_files)):
        if not input_files[x].endswith(".db"):
            gffutils.create_db(input_files[x], dbfn=input_files[x] + '.db', force=True,
                               disable_infer_genes=op_infer_genes, disable_infer_transcripts=op_infer_transcripts)
            input_files[x] += '.db'

    refdb = gffutils.FeatureDB(reference)
    # ref_genes = refdb.features_of_type('gene', order_by='start')
    reference_exons = list(refdb.features_of_type('exon'))

    asmdb, query_exons = [], []
    for x in input_files:
        asmdb.append(gffutils.FeatureDB(x))
        query_exons.append(list(asmdb[-1].features_of_type('exon')))


gene_name_id = {}
ref_junc = IntervalTree()
ref_intron = IntervalTree()


def parse_reference():
    global reference_exons
    global time2
    global ref_chrom_dict
    gene_strand_count = {}
    tran_strand_count = {}
    ref_chrom_dict = {}
    ref_gene_dict = {}
    ref_transcript_dict = {}
    time2 = time.time()
    cnt = 0
    number = 0
    intron_start = 0
    intron_end = 0
    prev_transcript_id = None
    for exon in tqdm(reference_exons):
        exon = Exon(exon)
        exon_interval = Interval(exon.start, exon.end + 1)
        if prev_transcript_id != exon.transcript_id:
            number = 1
            prev_transcript_id = exon.transcript_id
        else:
            number += 1
        # if cnt % 2 == 1:
        #     intron_end = exon.start + 1
        #     if prev_tran == exon.transcript_id:
        #         ref_junc[intron_start:intron_end] = ''
        #     if prev_strand == exon.strand:
        #         ref_intron[intron_start:intron_end] = ''
        #     print intron_start, intron_end
        #     prev_strand = exon.strand
        #     intron_start = exon.end
        #     prev_tran = exon.transcript_id
        #     cnt += 1
        #
        # else:
        #     if cnt > 0:
        #         intron_end = exon.start + 1
        #         if prev_tran == exon.transcript_id:
        #             ref_junc[intron_start:intron_end] = ''
        #         if prev_strand == exon.strand:
        #             ref_intron[intron_start:intron_end] = ''
        #         print intron_start, intron_end, 'else'
        #     prev_strand = exon.strand
        #     intron_start = exon.end
        #     prev_tran = exon.transcript_id
        #     cnt += 1
        if exon.attributes.get('gene_name'):
            gene_name = exon.attributes.get('gene_name')[0]
            gene_name_id[exon.gene_id] = gene_name
        if exon.chrom not in ref_chrom_dict:
            ref_chrom_dict[exon.chrom] = {}
        if exon.strand not in ref_chrom_dict[exon.chrom]:
            ref_intervals = IntervalTree()
            ref_exons_dict = {}
            ref_intervals[exon.start: exon.end + 1] = GtfInterval(exon.id, exon.transcript_id, exon.gene_id, number)
            ref_exons_dict[exon.id] = GtfExon(exon_interval, exon.transcript_id, exon.gene_id, number)

        else:
            ref_intervals, ref_exons_dict = ref_chrom_dict[exon.chrom][exon.strand]
            found = False
            dz = None
            for b, e, z in ref_intervals[exon_interval]:
                if b == exon_interval.begin and e == exon_interval.end:
                    found = True
                    dz = z
                    break
            if not found:
                ref_intervals[exon.start: exon.end + 1] = GtfInterval(exon.id, exon.transcript_id, exon.gene_id, number)
            else:
                if dz is None:
                    print 'Interval not found'
                    exit(1)
                dz.add(exon.id, exon.transcript_id, exon.gene_id, number)

            if exon.id not in ref_exons_dict:
                ref_exons_dict[exon.id] = GtfExon(exon_interval, exon.transcript_id, exon.gene_id, number)
            else:
                ref_exons_dict[exon.id].add(exon_interval, exon.transcript_id, exon.gene_id, number)

        temp = [ref_intervals, ref_exons_dict]
        ref_chrom_dict[exon.chrom][exon.strand] = temp

        if exon.transcript_id not in ref_transcript_dict:
            ref_transcript_dict[exon.transcript_id] = GtfTranscript(exon_interval, exon.id, exon.gene_id)
        else:
            ref_transcript_dict[exon.transcript_id].add(exon_interval, exon.id, exon.gene_id)

        if exon.gene_id not in ref_gene_dict:
            ref_gene_dict[exon.gene_id] = GtfGene(exon_interval, exon.id, exon.transcript_id)
        else:
            ref_gene_dict[exon.gene_id].add(exon_interval, exon.id, exon.transcript_id)

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

    for chromo in ref_chrom_dict:
        for strand in ref_chrom_dict[chromo]:
            ref_chrom_dict[chromo][strand].append(dict())
            ref_chrom_dict[chromo][strand].append(dict())

    for chromo in tran_strand_count:
        for tran, strands in tran_strand_count[chromo].items():
            if strands['+'] > strands['.'] and strands['+'] > strands['-']:
                strand = '+'
            elif strands['-'] > strands['.'] and strands['-'] > strands['+']:
                strand = '-'
            else:
                strand = '.'
            ref_chrom_dict[chromo][strand][2][tran] = ref_transcript_dict[tran]

    for chromo in gene_strand_count:
        for gene, strands in gene_strand_count[chromo].items():
            if strands['+'] > strands['.'] and strands['+'] > strands['-']:
                strand = '+'
            elif strands['-'] > strands['.'] and strands['-'] > strands['+']:
                strand = '-'
            else:
                strand = '.'
            ref_chrom_dict[chromo][strand][3][gene] = ref_gene_dict[gene]

    # return ref_chrom_dict


def parse_query(num):
    exonbest = open(input_files[num].split('/')[-1] + '.exon', 'w')
    exonbest.write("ExonID\tQuery(Cordinates[strand]|Transcript[exon_number])\tMatch_Type\t" +
                   "Reference(Chromosome:Best_Match_Cordinates|Transcript[exon_number])" +
                   "\tShared\tBase_Difference\n")
    gene_strand_count = {}
    tran_strand_count = {}
    quer_chrom_dict = {}
    quer_transcript_dict = {}
    quer_gene_dict = {}
    cnt = 0
    number = 0
    prev_transcript_id = None
    for exon in tqdm(query_exons[num]):
        cnt += 1
        exon = Exon(exon)
        exon_interval = Interval(exon.start, exon.end + 1)
        if prev_transcript_id != exon.transcript_id:
            number = 1
            prev_transcript_id = exon.transcript_id
        else:
            number += 1
        if exon.chrom not in quer_chrom_dict:
            quer_chrom_dict[exon.chrom] = {}

        if exon.strand not in quer_chrom_dict[exon.chrom]:
            quer_intervals = IntervalTree()
            quer_exons_dict = {}
            quer_intervals[exon.start: exon.end + 1] = GtfInterval(exon.id, exon.transcript_id, exon.gene_id, number)
            quer_exons_dict[exon.id] = GtfExon(exon_interval, exon.transcript_id, exon.gene_id, number)

        else:
            quer_intervals, quer_exons_dict = quer_chrom_dict[exon.chrom][exon.strand]
            found = False
            dz = None
            for b, e, z in quer_intervals[exon_interval]:
                if b == exon_interval.begin and e == exon_interval.end:
                    found = True
                    dz = z
                    break
            if not found:
                quer_intervals[exon.start: exon.end + 1] = GtfInterval(exon.id, exon.transcript_id, exon.gene_id, number)
            else:
                if dz is None:
                    print 'Interval not found'
                    exit(1)
                dz.add(exon.id, exon.transcript_id, exon.gene_id, number)

            if exon.id not in quer_exons_dict:
                quer_exons_dict[exon.id] = GtfExon(exon_interval, exon.transcript_id, exon.gene_id, number)
            else:
                quer_exons_dict[exon.id].add(exon_interval, exon.transcript_id, exon.gene_id, number)

        temp = [quer_intervals, quer_exons_dict]
        quer_chrom_dict[exon.chrom][exon.strand] = temp
        # bests = match_interval(next(iter(quer_intervals[exon.start: exon.end + 1])), exon.chrom, exon.strand)
        for cinter in quer_intervals[exon_interval]:
            if cinter.begin == exon_interval.begin and cinter.end == exon_interval.end:
                bests = match_interval(cinter, exon.chrom, exon.strand)
        if bests:
            for bintr, bval in bests.items():
                exonbest.write("Exon_" + str(cnt) + "\t" + str(exon.start) + "-" + str(exon.end) + "[" + exon.strand +
                               "]|" + str(exon.transcript_id) + "[" + str(number) + "]" + "\t" + bval[2] + "\t"
                               + exon.chrom + ":" + str(bintr.begin) + "-" + str(bintr.end - 1) + "[" + bval[4] + "]|"
                               + bintr.data.tran_number + "\t" + str(bval[0]) + "\t(" + str(bintr.begin - exon.start)
                               + "," + str(exon.end - bintr.end + 1) + ")\n")
        else:
            exonbest.write(
                "Exon_" + str(cnt) + "\t" + str(exon.start) + "-" + str(exon.end) + "[" + exon.strand + "]|" +
                str(exon.transcript_id) + "[" + str(number) + "]" + "\tNovel" + "\t-" * 3 + "\n")

        if exon.transcript_id not in quer_transcript_dict:
            quer_transcript_dict[exon.transcript_id] = GtfTranscript(exon_interval, exon.id, exon.gene_id)
        else:
            quer_transcript_dict[exon.transcript_id].add(exon_interval, exon.id, exon.gene_id)

        if exon.gene_id not in quer_gene_dict:
            quer_gene_dict[exon.gene_id] = GtfGene(exon_interval, exon.id, exon.transcript_id)
        else:
            quer_gene_dict[exon.gene_id].add(exon_interval, exon.id, exon.transcript_id)

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

    for chromo in quer_chrom_dict:
        for strand in quer_chrom_dict[chromo]:
            quer_chrom_dict[chromo][strand].append(dict())
            quer_chrom_dict[chromo][strand].append(dict())

    for chromo in tran_strand_count:
        for tran, strands in tran_strand_count[chromo].items():
            if strands['+'] > strands['.'] and strands['+'] > strands['-']:
                strand = '+'
            elif strands['-'] > strands['.'] and strands['-'] > strands['+']:
                strand = '-'
            else:
                strand = '.'
            quer_chrom_dict[chromo][strand][2][tran] = quer_transcript_dict[tran]

    for chromo in gene_strand_count:
        for gene, strands in gene_strand_count[chromo].items():
            if strands['+'] > strands['.'] and strands['+'] > strands['-']:
                strand = '+'
            elif strands['-'] > strands['.'] and strands['-'] > strands['+']:
                strand = '-'
            else:
                strand = '.'
            quer_chrom_dict[chromo][strand][3][gene] = quer_gene_dict[gene]

    return quer_chrom_dict


def get_all(chrom_dict, type_, interval):
    res = chrom_dict[0][interval]
    if len(res) > 0:
        temp = set()
        if type_ == 'transcripts':
            for intr in res:
                temp |= intr.data.transcriptIds
        elif type_ == 'genes':
            for intr in res:
                temp |= intr.data.geneIds
        return temp
    return None


def is_overlapped(feature_1, feature_2):
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
    end_1 -= 1
    end_2 -= 1
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
#     print overlaps, '\nafter\n'
#     overlaps = [v for x in overlaps for k, v in x.items() if v in ref_iso_intervals]
#     print overlaps
#     exit()


def find_super_loci(quer_chrom_dict):
    global ref_chrom_dict
    loc_ind = {}
    ind = 0
    linked = {}
    loci = {}
    super_loci = {}
    quer_loci_gene = {}
    for chrom in quer_chrom_dict:
        for strand in quer_chrom_dict[chrom]:
            quer_intervals, quer_exons_dict, quer_transcript_dict, quer_gene_dict =\
                quer_chrom_dict[chrom][strand]
            for locus in quer_gene_dict:
                quer = locus
                locus = quer_gene_dict[locus]
                interval = Interval(locus.start, locus.end - 1)
                found_overlap = False
                if chrom in ref_chrom_dict:
                    if strand in ref_chrom_dict[chrom]:
                    # for st in ref_chrom_dict[chrom]:
                        ref = get_all(ref_chrom_dict[chrom][strand], 'genes', interval)
                        if ref:
                            # if st != '.' and strand != '.':
                            #     if st != strand:
                            #         print 'found a locus overlap on opposite strand reference strand is ' + st \
                            #               + ' interval is ' + str(interval) + ' ' + strand + '\n', ref
                            found_overlap = True
                            for loc in ref:
                                if loc in loc_ind:
                                    prev_ind = loc_ind[loc]
                                    break
                                else:
                                    prev_ind = 0
                            if prev_ind == 0:
                                ind += 1
                                linked[ind] = [set(), set()]
                            for loc in ref:
                                if prev_ind != 0:
                                    loc_ind[loc] = prev_ind
                                    linked[prev_ind][1].add(loc)
                                    linked[prev_ind][0].add(quer)
                                else:
                                    loc_ind[loc] = ind
                                    linked[ind][1].add(loc)
                                    linked[ind][0].add(quer)

                            if quer not in quer_loci_gene:
                                quer_loci_gene[quer] = [chrom, strand]

                if not found_overlap:
                    if quer not in loci:
                        loci[quer] = [chrom, strand, interval]
                    else:
                        if loci[quer][0] != chrom or loci[quer][1] != strand or loci[quer][2] != interval:
                            print "Something went wrong dual gene", quer, interval
                            exit(1)
    for v in linked.values():
        v[0] = list(v[0])
        if v[0]:
            # print v, ' v\nv0 ', v[0], '\nv1 ', v[1]
            super_loci[tuple(v[0])] = quer_loci_gene[v[0][0]] + list(v[1])
    return [loci, super_loci]


def write_loci_file(quer_chrom_dict, loci, super_loci, x):
    formatted = {}
    super_loci_counter = 1
    with open('gtfcomp.' + input_files[x].split('/')[-1] + '.loci', 'w') as f:
        for k, v in sorted(loci.items()):
            chrom, strand, interval = v
            line = k+'\t'+chrom+'['+strand+']'+str(interval.begin)+'-'+str(interval.end)+'\t-\t'
            # f.write(line)
            quer_gene_dict = quer_chrom_dict[chrom][strand][3]
            quer = sorted(quer_gene_dict[k].transcriptIds, reverse=True)
            for l in range(len(quer) - 1):
                line += quer.pop() + ','
                # f.write(quer.pop() + ',')
            # f.write(quer.pop() + '\n')
            line += quer.pop() + '\n'
            if chrom not in formatted:
                formatted[chrom] = {}
            if strand not in formatted[chrom]:
                formatted[chrom][strand] = {}
            if interval.begin not in formatted[chrom][strand]:
                formatted[chrom][strand][interval.begin] = [line]
            else:
                formatted[chrom][strand][interval.begin].append(line)
        for keys, v in super_loci.items():
            quer = set()
            chrom, strand, genes = v[0], v[1], v[2:]
            quer_gene_dict = quer_chrom_dict[chrom][strand][3]
            ref = ''
            start = 1000**1000
            end = 0
            for k in keys:
                if k in quer_gene_dict:
                    gen = quer_gene_dict[k]
                    quer |= gen.transcriptIds
                    start = min(start, gen.start)
                    end = max(end, gen.end)
            if quer:
                quer = sorted(quer, reverse=True)
                k = 'SLOC_'+str(super_loci_counter).zfill(6)
                super_loci_counter += 1
                ref_gene_dict = ref_chrom_dict[chrom][strand][3]
                for g in genes:
                    if g in ref_gene_dict:
                        gen = ref_gene_dict[g]
                        g = gene_name_id.get(g, g)
                        start = min(start, gen.start)
                        end = max(end, gen.end)
                        ref += g + '|'
                        trans = sorted(gen.transcriptIds, reverse=True)
                        for l in range(len(trans) - 1):
                            ref += trans.pop() + '|'
                        ref += trans.pop() + ','
                if ref:
                    ref = ref[:-1]
                    line = k+'\t'+chrom+'['+strand+']'+str(start)+'-'+str(end - 1)+'\t'+ref+'\t'
                    # f.write(line)
                    for l in range(len(quer) - 1):
                        line += quer.pop() + ','
                        # f.write(quer.pop() + ',')
                    # f.write(quer.pop() + '\n')
                    line += quer.pop() + '\n'
                    if chrom not in formatted:
                        formatted[chrom] = {}
                    if strand not in formatted[chrom]:
                        formatted[chrom][strand] = {}
                    if start not in formatted[chrom][strand]:
                        formatted[chrom][strand][start] = [line]
                    else:
                        formatted[chrom][strand][start].append(line)
        for ch in sorted(formatted):
            for st in sorted(formatted[ch]):
                for intr in sorted(formatted[ch][st]):
                    for item in formatted[ch][st][intr]:
                        f.write(item)


interval_matches = {}
interval_best_matches = {}
ref_matches = {}
ref_best_matches = {}


def match_interval(interval, chrom, stran):
    global interval_matches
    global interval_best_matches
    global ref_matches
    global ref_best_matches
    global ref_chrom_dict

    if interval in interval_matches:
        return interval_best_matches[interval]
    found = False
    if chrom in ref_chrom_dict:
        if stran in ref_chrom_dict[chrom]:
            found = ref_chrom_dict[chrom][stran][0][interval]
            ref_strand = stran
        if not found and stran != '.':
            if '.' in ref_chrom_dict[chrom]:
                found = ref_chrom_dict[chrom]['.'][0][interval]
                ref_strand = '.'
            elif stran == '+' and '-' in ref_chrom_dict[chrom]:
                found = ref_chrom_dict[chrom]['-'][0][interval]
                ref_strand = '-'
                # if found:
                    # print 'found match for ', interval, ' on reverse strand on reference'
            elif stran == '-' and '+' in ref_chrom_dict[chrom]:
                found = ref_chrom_dict[chrom]['+'][0][interval]
                ref_strand = '+'
                # if found:
                    # print 'found match for ', interval, ' on reverse strand on reference'
        elif not found:
            ref_strand = '.'
            if '+' in ref_chrom_dict[chrom]:
                found_pos = ref_chrom_dict[chrom]['+'][0][interval]
            if '-' in ref_chrom_dict[chrom]:
                found_neg = ref_chrom_dict[chrom]['-'][0][interval]
            if '.' in ref_chrom_dict[chrom]:
                found_un = ref_chrom_dict[chrom]['.'][0][interval]

            if not found_pos and not found_un and found_neg:
                found = found_neg
                ref_strand = '-'
            elif not found_neg and not found_un and found_pos:
                found = found_pos
                ref_strand = '+'
            elif not found_neg and not found_pos and found_un:
                found = found_pos
            else:
                found = found_pos | found_neg | found_un
                # comp = [x for x in found_pos]
                # pos_max = map(is_overlapped, comp, itertools.repeat(interval, len(comp)))
                # ii, pos_max = max(enumerate(pos_max), key=operator.itemgetter(1))
                # p_size = comp[ii].end - comp[ii].begin
                #
                # comp = [x for x in found_neg]
                # neg_max = map(is_overlapped, comp, itertools.repeat(interval, len(comp)))
                # ii, neg_max = max(enumerate(neg_max), key=operator.itemgetter(1))
                # n_size = comp[ii].end - comp[ii].begin
                # if (pos_max > neg_max) or (pos_max == neg_max and p_size <= n_size):
                #     found = found_pos
                # else:
                #     found = found_neg
                #
        if found:
            if interval.end == 126343416:
                print "GOT FOUND"
            comp = [x for x in found]
            comp_overlaps = map(is_overlapped, comp, itertools.repeat(interval, len(comp)))
            comp_judgements = [x[1] for x in comp_overlaps]
            comp_overlaps = [x[0] for x in comp_overlaps]
            comp_sizes = map(lambda xx: xx.end - 1 - xx.begin, comp)
            # comp = [Interval(intr.begin, intr.end) for intr in comp]
            # ii, ovlp = max(enumerate(comp_overlaps), key=operator.itemgetter(1))
            ii, ovlp, ml = [0], 0, False
            for i, v in enumerate(comp_overlaps):
                if v > ovlp or (v == ovlp and comp_sizes[i] < comp_sizes[ii[0]]):
                    ovlp = v
                    ii = [i]
                elif v == ovlp and comp_sizes[i] == comp_sizes[ii[0]]:
                    ii.append(i)
            # stinterval = str(interval.begin) + ':' + str(interval.end - 1)
            intr_size = interval.end - 1 - interval.begin

            if interval not in interval_best_matches:
                interval_best_matches[interval] = dict()
                for xi in ii:
                    interval_best_matches[interval][comp[xi]] = (ovlp, comp_sizes[xi],
                                                                 comp_judgements[xi], stran, ref_strand)
            else:
                print 'Dual Interval found ', interval, 'more than once ', interval_best_matches[interval],\
                    '\n', (comp[ii[0]], ovlp, comp_sizes[ii[0]], comp_judgements[ii[0]], stran, ref_strand)

            for i, v in enumerate(comp):
                if v in ref_best_matches:
                    ri = ref_best_matches[v]
                    rii = ri.itervalues().next()
                    if (comp_overlaps[i] > rii[0]) or (comp_overlaps[i] == rii[0] and intr_size < rii[1]):
                        ri = dict()
                        ri[interval] = (comp_overlaps[i], intr_size, comp_judgements[i], stran, ref_strand)
                    elif comp_overlaps[i] == rii[0] and comp_sizes[i] == rii[1]:
                        ri[interval] = (comp_overlaps[i], intr_size, comp_judgements[i], stran, ref_strand)
                else:
                    ref_best_matches[v] = dict()
                    ri = ref_best_matches[v]
                    ri[interval] = (comp_overlaps[i], intr_size, comp_judgements[i], stran, ref_strand)
                if v not in ref_matches:
                    ref_matches[v] = {}
                ref_matches[v][interval] = (comp_overlaps[i], intr_size, comp_judgements[i], stran, ref_strand)
                if interval not in interval_matches:
                    interval_matches[interval] = {}
                interval_matches[interval][v] = (comp_overlaps[i], comp_sizes[i],
                                                 comp_judgements[i], stran, ref_strand)
            return interval_best_matches[interval]


argument_parser()
parse_database()
parse_reference()
# quer = dict()
# quer['chr1'] = {}
# quer['chr1']['.'] = [IntervalTree()]
# quer['chr1']['.'][0][29560:29565] = 'asd'
#
# match_exons(quer)
for x in range(len(input_files)):
    quer_chrom_dict = parse_query(x)
    exonbest = open(input_files[x].split('/')[-1] + '.etracking', 'w')
    exonbest.write("ID\tQuery(Cordinates[strand]|Transcript[exon_number])\t"
                   "Reference(Ref_Cordinates[strand]|Transcript[exon_number])\t")
    exonbest.write("Type\t(Start,End)\tQuery_Stats\tReference_Stats\n")
    cnt = 0
    reported = set()
    for chrom in quer_chrom_dict:
        for strand in quer_chrom_dict[chrom]:
            for eintr in sorted(quer_chrom_dict[chrom][strand][0]):
                cnt += 1
                if eintr in reported:
                    continue
                elif eintr in interval_matches:
                    reported |= set([eintr])
                    qset = dict()
                    rset = dict()
                    done = False

                    for match in interval_matches[eintr]:
                        rset[match] = interval_matches[eintr][match]
                        for qmatch in ref_matches[match]:
                            qset[qmatch] = ref_matches[match][qmatch]
                    temp_qset = qset
                    while not done:
                        temp_rdict = dict()
                        temp_qdict = dict()
                        done = True
                        for q in temp_qset:
                            for match in interval_matches[q]:
                                if match not in rset:
                                    temp_rdict[match] = interval_matches[q][match]
                                    for qmatch in ref_matches[match]:
                                        if qmatch not in qset:
                                            done = False
                                            temp_qdict[qmatch] = ref_matches[match][qmatch]
                        temp_qset = temp_qdict
                        if temp_qdict:
                            qset.update(temp_qdict)
                        if temp_rdict:
                            rset.update(temp_rdict)

                    qset = sorted(qset)
                    rset = sorted(rset)

                    if len(qset) == 1 and len(rset) == 1:
                        mtype = "SQSR"
                        rm = next(iter(rset))
                        qm = next(iter(qset))
                        cord = "{}-{}[{}]|{}\t".format(eintr.begin, eintr.end - 1, strand,
                                                            "|".join(sorted(eintr.data.transcriptIds)))
                        cord += "r|{}-{}[{}]|{}\t".format(rm.begin, rm.end - 1, ref_matches[rm][qm][4],
                                                            "|".join(sorted(rm.data.transcriptIds)))
                        stats = "({},{})\t-\t-".format(rm.begin - eintr.begin, eintr.end - rm.end)

                    elif len(qset) == 1 and len(rset) > 1:
                        mtype = "SQMR"
                        start = 10000*10000
                        end = 0
                        qm = next(iter(qset))
                        cord = "{}-{}[{}]|{}\t".format(eintr.begin, eintr.end - 1, strand,
                                                       "|".join(sorted(eintr.data.transcriptIds)))
                        stat = ""
                        began = False
                        for bintr in rset:
                            cord += "r|{}-{}[{}]|{}, ".format(bintr.begin, bintr.end - 1, ref_matches[bintr][qm][4],
                                                              "|".join(sorted(bintr.data.transcriptIds)))
                            if began:
                                stat += str(prev_end - bintr.begin) + " "
                            else:
                                began = True
                            start = min(start, bintr.begin)
                            end = max(end, bintr.end)
                            prev_end = bintr.end - 1
                        cord = cord[:-2]
                        stats = "({},{})\t-\t[{}]".format(start - eintr.begin, eintr.end - end, stat)

                    elif len(qset) > 1 and len(rset) == 1:
                        mtype = "MQSR"
                        start = 10000 * 10000
                        end = 0
                        rm = next(iter(rset))
                        cord = ""
                        stat = ""
                        began = False
                        for bintr in qset:
                            reported |= set([bintr])
                            cord += "{}-{}[{}]|{}, ".format(bintr.begin, bintr.end - 1, interval_matches[bintr][rm][3],
                                                            "|".join(sorted(bintr.data.transcriptIds)))
                            if began:
                                stat += str(prev_end - bintr.begin) + " "
                            else:
                                began = True
                            start = min(start, bintr.begin)
                            end = max(end, bintr.end)
                            prev_end = bintr.end - 1
                        cord = cord[:-2]
                        cord += "\tr|{}-{}[{}]|{}".format(rm.begin, rm.end - 1, ref_matches[rm][qset[0]][4],
                                                       "|".join(sorted(rm.data.transcriptIds)))
                        stats = "({},{})\t[{}]\t-".format(rm.begin - start, end - rm.end, stat)

                    elif len(qset) > 1 and len(rset) > 1:
                        mtype = "MQMR"
                        qstart = 10000 * 10000
                        qend = 0
                        cord = ""
                        qstat = ""
                        began = False
                        rgabs = {}
                        qgabs = {}
                        for bintr in qset:
                            reported.add(bintr)

                            cord += "{}-{}[{}]|{}, ".format(bintr.begin, bintr.end - 1, strand,
                                                            "|".join(sorted(bintr.data.transcriptIds)))
                            if began:
                                qstat += str(prev_end - bintr.begin) + " "
                                qgabs[str(bintr.begin)+":"+str(prev_end)] = str(bintr.begin - prev_end)
                            else:
                                began = True
                            qstart = min(qstart, bintr.begin)
                            qend = max(qend, bintr.end)
                            prev_end = bintr.end - 1
                        cord = cord[:-2]
                        rstat = ""
                        rstart = 10000 * 10000
                        rend = 0
                        began = False
                        for bintr in rset:
                            cord += "\tr|{}-{}[{}]|{}, ".format(bintr.begin, bintr.end - 1, strand,
                                                            "|".join(sorted(bintr.data.transcriptIds)))
                            if began:
                                rstat += str(prev_end - bintr.begin) + " "
                                rgabs[str(bintr.begin)+":"+str(prev_end)] = str(bintr.begin - prev_end)
                            else:
                                began = True
                            rstart = min(rstart, bintr.begin)
                            rend = max(rend, bintr.end)
                            prev_end = bintr.end - 1
                        cord = cord[:-2]
                        rstat = [rgabs[k] for k in sorted(rgabs.keys())]
                        qstat = [qgabs[k] for k in sorted(qgabs.keys())]
                        stats = "({},{})\t[{}]\t[{}]".format(rstart - qstart, qend - rend, " ".join(qstat),
                                                             " ".join(rstat))

                    exonbest.write("id_{}\t{}\t".format(cnt, cord))
                    exonbest.write("{}\t{}\n".format(mtype, stats))
                else:
                    exonbest.write("id_{}\t{}-{}[{}]|{}\t-\t".format(cnt, eintr.begin, eintr.end - 1, strand,
                                                                     ",".join(eintr.data.transcriptIds)))
                    exonbest.write("SQNR\t-\t-\t-\n")

    # loci, super_loci = find_super_loci(quer_chrom_dict)
    # write_loci_file(quer_chrom_dict, loci, super_loci, x)
print "total time ", time.time() - time1
print "time without gffutils ", time.time() - time2
