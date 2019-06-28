import itertools
import six

class MatchingDicts:
    interval_matches = {}
    interval_best_matches = {}
    ref_matches = {}
    ref_best_matches = {}

    def __init__(self):
        MatchingDicts.get_dicts()

    @staticmethod
    def get_dicts():
        return [MatchingDicts.interval_matches, MatchingDicts.interval_best_matches, MatchingDicts.ref_matches, MatchingDicts.ref_best_matches]

    @staticmethod
    def set_dicts(dicts):
        MatchingDicts.interval_matches, MatchingDicts.interval_best_matches, MatchingDicts.ref_matches, MatchingDicts.ref_best_matches = dicts

    @staticmethod
    def get_interval_matches():
        return MatchingDicts.interval_matches

    @staticmethod
    def get_interval_best_matches():
        return MatchingDicts.interval_best_matches

    @staticmethod
    def get_ref_matches():
        return MatchingDicts.ref_matches

    @staticmethod
    def get_ref_best_matches():
        return MatchingDicts.ref_best_matches

def is_overlapped(feature_1, feature_2):
    """Calculates the size of overlap and judgement of two features if exists"""

    judgement = 'No_Overlap'
    start_1 = feature_1.begin
    end_1 = feature_1.end

    start_2 = feature_2.begin
    end_2 = feature_2.end

    size_1 = end_1 - start_1 + 1
    size_2 = end_2 - start_2 + 1

    if start_1 > start_2 or (start_1 == start_2 and size_2 > size_1):
        start_1, end_1, size_1, start_2, end_2, size_2 = start_2, end_2, size_2, start_1, end_1, size_1

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
        return [end_1 - start_2 + 1, judgement]

    return [size_2, judgement]

def match_interval(quer_chrom_dict, ref_chrom_dict):
    """Matches query interval to the reference"""

    matching_dicts = MatchingDicts()

    interval_matches, interval_best_matches, ref_matches, ref_best_matches = matching_dicts.get_dicts()

    # Indicates if any match was found
    found = False

    # No notes initially
    note = 0
    for chrom in quer_chrom_dict:
        for stran in quer_chrom_dict[chrom]:
            for interval in quer_chrom_dict[chrom][stran][0]:
                if chrom in ref_chrom_dict:
                    if stran in ref_chrom_dict[chrom]:
                        # Search Same strand first
                        found = ref_chrom_dict[chrom][stran][0][interval.begin:interval.end]

                    # If nothing was found search unknown first then opposite strand (strand is either + or -)
                    if not found and stran != '.':
                        # Search unknown first
                        if '.' in ref_chrom_dict[chrom]:
                            found = ref_chrom_dict[chrom]['.'][0][interval.begin:interval.end]
                            # Note 1: A match was found on the unknown strand of the reference
                            note = 1
                        # Then search opposite
                        if not found:
                            if stran == '+' and '-' in ref_chrom_dict[chrom]:
                                found = ref_chrom_dict[chrom]['-'][0][interval.begin:interval.end]

                            elif stran == '-' and '+' in ref_chrom_dict[chrom]:
                                found = ref_chrom_dict[chrom]['+'][0][interval.begin:interval.end]

                            if found:
                                # Note 2: A match was found on the opposite strand on the reference
                                note = 2
                                # print('found match for ', interval, ' on reverse strand on reference')

                    # If strand was unknown search both strands
                    elif not found and stran == '.':
                        found_pos = set()
                        found_neg = set()

                        if '+' in ref_chrom_dict[chrom]:
                            found_pos = ref_chrom_dict[chrom]['+'][0][interval.begin:interval.end]

                        if '-' in ref_chrom_dict[chrom]:
                            found_neg = ref_chrom_dict[chrom]['-'][0][interval.begin:interval.end]

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
    matching_dicts.set_dicts([interval_matches, interval_best_matches, ref_matches, ref_best_matches])
    return matching_dicts.get_dicts()