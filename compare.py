import gffutils
import datetime

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
        self.__set_exon_cat()

    def __set_exon_cat(self):
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
                if longest_first_exon_id == exon_.id:
                    exon_.type_ = 'Type 1'
                else:
                    exon_.type_ = 'Type 1/2'
            elif exon_.cat_ == 'last':
                if longest_last_exon_id == exon_.id:
                    exon_.type_ = 'Type 2'
                else:
                    exon_.type_ = 'Type 2/2'

        for cexon in listofexons:
            if cexon.cat_ == 'mid':
                for nexon in listofexons:
                    if cexon.id != nexon.id:
                        if nexon.cat_ == 'first':
                            if Exon.match_exons(cexon, nexon)[0]:
                                cexon.type_ = 'Type 4'
                                break
                        elif nexon.cat_ == 'last':
                            if Exon.match_exons(cexon, nexon)[0]:
                                cexon.type_ = 'Type 5'
                                break
                if cexon.type_ is None:
                    cexon.type_ = 'Type 3'

    @staticmethod
    def match_exons(exon1, exon2):
        if exon1.start > exon2.start:
            temp = exon1
            exon1 = exon2
            exon2 = temp
        if exon1.end - exon2.start > 0:
            return [True, exon1.end - exon2.start]
        return [False, None]


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
    count, start = 0, 0
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
            overlap, basepair = Exon.match_exons(asm_feature, ref_feature)
            if overlap and basepair > max_basepair:
                matched_ref_feature = Feature(ref_feature)
                max_basepair = basepair
            elif overlap and basepair == max_basepair and max_basepair != 0:
                if ref_feature.size < matched_ref_feature.size:
                    matched_ref_feature = Feature(ref_feature)
            if overlap:
                if matched_ref_feature.id not in ref_matches.keys():
                    ref_matches[matched_ref_feature.id] = basepair
                else:
                    ref_matches[matched_ref_feature.id] = max(ref_matches[matched_ref_feature.id], basepair)
                if asm_feature.id not in asm_matches.keys():
                    asm_matches[asm_feature.id] = basepair
                else:
                    asm_matches[asm_feature.id] = max(asm_matches[asm_feature.id], basepair)

        if matched_ref_feature is not None:
            if matched_ref_feature.id not in best_matches.keys():
                best_matches[matched_ref_feature.id] = [(asm_feature.id, basepair)]
            else:
                best_matches[matched_ref_feature.id].append((asm_feature.id, basepair))
        # count += 1
        # print count
    return best_matches


refdb = gffutils.FeatureDB("databases/poly_brain_stringtie_merged.gtf.db")
ref_exons = list(refdb.features_of_type('exon', order_by='start'))


asmdb = gffutils.FeatureDB("databases/poly_brain_scallop_merged.gtf.db")
asm_exons = list(asmdb.features_of_type('exon', order_by='start'))

ref_genes = refdb.features_of_type('gene', order_by='start')
ref_parents = Parents(refdb)

d1 = datetime.datetime.now()

best_match = find_best_match(asm_exons, ref_exons)

d2 = datetime.datetime.now()

print (d2 - d1)

asm_parents = Parents(asmdb)

for gene in ref_genes:
    refexons = []
    child_exons = refdb.children(gene.id, featuretype='exon')
    for exon in child_exons:
        refexons.append(Exon(exon, ref_parents.transcript(exon)))
    Exon.set_exon_type(refexons)
    with open("ref_exon_table.gtf", 'a') as file_:
        for exon in refexons:
            file_.write("---- {} \t {} \t {} \t {} \t {} \t {} \t {} \t {}\n"
                        .format(ref_parents.gene(exon).id, exon.transcript.id, exon.id, exon.start,
                                exon.end, exon.size, exon.type_))
            if exon.id in best_match.keys():
                for match, bp in best_match[exon.id]:
                    asm_exon = Exon(asmdb[asmdb[match]], asm_parents.transcript(match))
                    if bp == asm_exon.size:
                        relation = 'Perfect'
                    else:
                        relation = 'Partial'
                    file_.write("**** {} \t {} \t {} \t {} \t {} \t {} \t {} \t {}\n"
                                .format(asm_parents.gene(asm_exon).id, asm_parents.transcript(asm_exon).id, asm_exon.id,
                                        asm_exon.start, asm_exon.end, asm_exon.size, asm_exon.cat_,    asm_exon.type_))
