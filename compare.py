import gffutils
from argparse import ArgumentParser


class Feature(gffutils.Feature):

    def __init__(self, feature):
        gffutils.Feature.__init__(self, feature.seqid, feature.source, feature.featuretype,
                                  feature.start, feature.end, feature.score, feature.strand, feature.frame,
                                  feature.attributes, feature.extra, feature.bin, feature.id, feature.dialect,
                                  feature.file_order, feature.keep_order, feature.sort_attribute_values)
        self.size = feature.end - feature.start + 1


class Exon(Feature):

    type_ = None
    cat_ = None

    def __init__(self, feature, transcript):
        Feature.__init__(self, feature)
        self.transcript = transcript

    def set_exon_cat(self):
        if self.start == self.transcript.start:
            self.cat_ = 'first'
        elif self.end == self.transcript.end:
            self.cat_ = 'last'
        else:
            self.cat_ = 'mid'

    @staticmethod
    def set_exon_type(listofexons):
        max_first_exon = 0
        max_last_exon = 0
        longest_first_exon_id = None
        longest_last_exon_id = None

        for exon_ in listofexons:
            if exon_.cat_ == 'first':
                if exon_.size > max_first_exon:
                    max_first_exon = exon_.size
                    longest_first_exon_id = exon_.id
            elif exon_.cat_ == 'last':
                if exon_.size > max_last_exon:
                    max_last_exon = exon_.size
                    longest_last_exon_id = exon_.id

        for exon_ in listofexons:
            if exon_.cat_ == 'first':
                if exon_.id == longest_first_exon_id:
                    exon_.type_ = 'Type 1'
                else:
                    exon_.type_ = 'Type 1/2'
            elif exon_.cat_ == 'last':
                if exon_.id == longest_last_exon_id:
                    exon_.type_ = 'Type 2'
                else:
                    exon_.type_ = 'Type 2/2'

        for cexon in listofexons:
            if cexon.cat_ == 'mid':
                for nexon in listofexons:
                    if cexon.id != nexon.id:
                        if nexon.cat_ == 'first':
                            basepair = Exon.match_exons(cexon, nexon)
                            if basepair > threshold:
                                cexon.type_ = 'Type 4'
                                nexon.type_ = 'Type 4'
                                break
                        elif nexon.cat_ == 'last':
                            basepair = Exon.match_exons(cexon, nexon)
                            if basepair > threshold:
                                cexon.type_ = 'Type 5'
                                nexon.type_ = 'Type 5'
                                break
                if cexon.type_ is None:
                    cexon.type_ = 'Type 3'

    @staticmethod
    def match_exons(exon1, exon2):
        if exon1.start > exon2.start or (exon1.start == exon2.start and exon2.size > exon1.size):
            temp_ = exon1
            exon1 = exon2
            exon2 = temp_
        if exon2.end > exon2.end:
            return exon1.end - exon2.start
        return exon2.size


class Parents:
    def __init__(self, db):
        self.db = db

    def gene(self, feature):
        return next(self.db.parents(feature.id, featuretype='gene'))

    def transcript(self, feature):
        return next(self.db.parents(feature.id, featuretype='transcript'))

    def transcript(self, feature):
        return next(self.db.parents(feature, featuretype='transcript'))

    def transcriptstart(self, feature):
        return next(self.db.parents(feature.id, featuretype='transcript')).start

    def transcriptend(self, feature):
        return next(self.db.parents(feature.id, featuretype='transcript')).end


def find_best_match(asm_features, ref_features):
    best_matches = {}
    ref_matches = {}
    asm_matches = {}
    start = 0
    sz = len(asm_features)
    for index in range(sz):
        asm_feature = Feature(asm_features[index])
        matched_ref_feature, max_basepair = None, 0
        for refind in range(start, len(ref_features)):
            ref_feature = Feature(ref_features[refind])
            if asm_feature.end < ref_feature.start:
                break
            if asm_feature.start > ref_feature.end:
                continue
            if index != sz - 1:
                if ref_feature.end < asm_features[index + 1].start:
                    start = refind
            basepair = Exon.match_exons(asm_feature, ref_feature)
            if basepair > max_basepair:
                matched_ref_feature = Feature(ref_feature)
                max_basepair = basepair
            elif basepair == max_basepair and max_basepair != 0:
                if ref_feature.size < matched_ref_feature.size:
                    matched_ref_feature = Feature(ref_feature)
            if basepair > 0:
                if ref_feature.id not in ref_matches.keys():
                    ref_matches[ref_feature.id] = basepair
                else:
                    ref_matches[ref_feature.id] = max(ref_matches[ref_feature.id], basepair)
                if asm_feature.id not in asm_matches.keys():
                    asm_matches[asm_feature.id] = basepair
                else:
                    asm_matches[asm_feature.id] = max(asm_matches[asm_feature.id], basepair)

        if matched_ref_feature is not None:
            if matched_ref_feature.id not in best_matches.keys():
                best_matches[matched_ref_feature.id] = [(asm_feature.id, max_basepair)]
            else:
                best_matches[matched_ref_feature.id].append((asm_feature.id, max_basepair))
    return [best_matches, ref_matches, asm_matches]


def write_junctions(ref=None, asm=None):
    if ref is not None:
        junctions = str(ref.start)
        if ref.type_ != "Type 2":
            junctions += "\t" + ref.type_
        else:
            junctions += "\t" + "Type 3"
        if asm is not None:
            if ref.start == asm.start:
                relation_ = "Precise"
            else:
                relation_ = "Imprecise"
            junctions += "\t;\t" + str(asm.start) + "\t" + relation_ + "\n"
        else:
            junctions += "\tNot Found\t" + "\n"
    else:
        junctions = str(asm.start) + "\t Novel\n"

    if ref is not None:
        junctions += str(ref.end)
        if ref.type_ != "Type 1":
            junctions += "\t" + ref.type_
        else:
            junctions += "\t" + "Type 3"
        if asm is not None:
            if ref.end == asm.end:
                relation_ = "Precise"
            else:
                relation_ = "Imprecise"
            junctions += "\t;\t" + str(asm.end) + "\t" + relation_ + "\n"
        else:
            junctions += "\tNot Found" + "\n"
    else:
        junctions += str(asm.end) + "\t Novel\n"

    return junctions


def argument_parser():
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
    return [input_annotations, args.refgtf, args.th[0], args.it, args.ig]


arguments = argument_parser()
input_files = arguments[0]
reference = str(arguments[1])
threshold = arguments[2]
op_infer_transcripts = arguments[3]
op_infer_genes = arguments[4]
del arguments

if not reference.endswith(".db"):
    gffutils.create_db(reference, dbfn=reference+'.db', force=True,
                       disable_infer_genes=op_infer_genes, disable_infer_transcripts=op_infer_transcripts)
    reference += '.db'

for x in range(len(input_files)):
    if not input_files[x].endswith(".db"):
        gffutils.create_db(input_files[x], dbfn=input_files[x]+'.db', force=True,
                           disable_infer_genes=op_infer_genes, disable_infer_transcripts=op_infer_transcripts)
        input_files[x] += '.db'

refdb = gffutils.FeatureDB(reference)
ref_genes = refdb.features_of_type('gene', order_by='start')
ref_exons = list(refdb.features_of_type('exon', order_by='start'))
ref_parents = Parents(refdb)

asmdb, asm_exons, asm_parents = [], [], []
best_match, refmatches, asmmatches = [], [], []
file_table, file_junction = [], []


for x in input_files:
    asmdb.append(gffutils.FeatureDB(x))
    asm_exons.append(list(asmdb[-1].features_of_type('exon', order_by='start')))
    asm_parents.append(Parents(asmdb[-1]))
    temp = find_best_match(asm_exons[-1], ref_exons)
    best_match.append(temp[0])
    refmatches.append(temp[1])
    asmmatches.append(temp[2])
    del temp

    file_table.append(open(x + ".comp", 'w'))
    file_junction.append(open(x + ".junc", 'w'))

for gene in ref_genes:
    refexons = []
    child_exons = refdb.children(gene.id, featuretype='exon')
    for exon in child_exons:
        refexons.append(Exon(exon, ref_parents.transcript(exon)))
        refexons[-1].set_exon_cat()
    Exon.set_exon_type(refexons)
    for exon in refexons:
        for ii in range(len(input_files)):
            file_table[ii].write("-ref {} \t {} \t {} \t {} \t {} \t {} \t {}"
                                 .format(ref_parents.gene(exon).id, exon.transcript.id, exon.id, exon.start,
                                         exon.end, exon.size, exon.type_))
            if exon.id in refmatches[ii].keys():
                if exon.id in best_match[ii].keys():
                    for match, bp in best_match[ii][exon.id]:
                        matched = False
                        asm_exon = Exon(asmdb[ii][match], asm_parents[ii].transcript(match))
                        if bp > threshold:
                            matched = True
                            if bp == (asm_exon.size + exon.size)/2:
                                relation = 'Perfect'
                            else:
                                relation = 'Partial'
                            file_table[ii].write("\n*asm {} \t {} \t {} \t {} \t {} \t {} \t {}"
                                             .format(asm_parents[ii].gene(asm_exon).id, asm_exon.transcript.id,
                                                     asm_exon.id, asm_exon.start, asm_exon.end, asm_exon.size, relation))
                            file_junction[ii].write(write_junctions(exon, asm_exon))
                    if not matched:
                        file_table[ii].write("\t Not Found I")
                        file_junction[ii].write(write_junctions(ref=exon))
                else:
                    bp = refmatches[ii][exon.id]
                    if bp > threshold:
                        file_table[ii].write("\t Not Found III")
                    else:
                        file_table[ii].write("\t Not Found II")
                    file_junction[ii].write(write_junctions(ref=exon))
            else:
                file_table[ii].write("\t Not Found I")
                file_junction[ii].write(write_junctions(ref=exon))
            file_table[ii].write("\n")

for x in range(len(input_files)):
    file_table[ii].write("\n----------------------------Novel Assembled Exons----------------------------\n")
    file_junction[ii].write("\n----------------------------Novel Assembled Junctions----------------------------\n")
    for exon in asm_exons[ii]:
        exon = Exon(exon, asm_parents[ii].transcript(exon))
        if exon.id not in asmmatches[ii].keys():
            file_table[ii].write("**** {} \t {} \t {} \t {} \t {} \t {} \t {}\n"
                             .format(asm_parents[ii].gene(exon).id, exon.transcript.id, exon.id,
                                     exon.start, exon.end, exon.size, "Novel I"))
            file_junction[ii].write(write_junctions(asm=exon))
        else:
            bp = asmmatches[ii][exon.id]
            if bp < threshold:
                file_table[ii].write("**** {} \t {} \t {} \t {} \t {} \t {} \t {}\n"
                                 .format(asm_parents[ii].gene(exon).id, exon.transcript.id, exon.id,
                                         exon.start, exon.end, exon.size, "Novel II"))
                file_junction[ii].write(write_junctions(asm=exon))
                
for x in range(len(input_files)):
    file_table[ii].close()
    file_junction[ii].close()
    

