from interval_matching import MatchingDicts
from notes import NOTES
import six

def write_exon(input_file, gtf_exons):
    """Writes the .exon file"""

    dotexon = open(input_file + '.exon', 'w')
    if input_file == 'ref':
        interval_best_matches = MatchingDicts.ref_best_matches
        dotexon.write("ExonID\tChromosome\tReference(Coordinates[strand]|Transcript[exon_number])\tMatch_Type\t" +
                "Query(Best_Match_Coordinates|Transcript[exon_number])\tShared\tBase_Difference\tNotes\n")
    else:
        interval_best_matches = MatchingDicts.interval_best_matches
        dotexon.write("ExonID\tChromosome\tQuery(Coordinates[strand]|Transcript[exon_number])\tMatch_Type\t" +
                    "Reference(Best_Match_Coordinates|Transcript[exon_number])\tShared\tBase_Difference\tNotes\n")

    for cnt, exon_id in enumerate(gtf_exons):
        exon = gtf_exons[exon_id]
        interval_begin, interval_end = exon.begin, exon.end
        cinter = exon.gtf_interval
        bests = interval_best_matches[interval_begin, interval_end]
        # If a match (best match) was found write each match in .exon file
        if bests:
            for bintr, bval in bests.items():
                dotexon.write('{}\t{}\t{}-{}[{}]|{}[{}]\t{}\t{}-{}[{}]|{}\t{}\t({},{})\t({})\n'.format(
                    cnt + 1, exon.chrom, cinter.begin, cinter.end - 1, cinter.data.strand, exon.transcript_id,
                    cinter.data.transcriptIds[exon.transcript_id], bval[1], bintr.begin, bintr.end - 1,
                    bintr.data.strand, '|'.join(['{}[{}]'.format(k, v) for k, v in bintr.data.transcriptIds.items()]),
                    bval[0], bintr.begin - cinter.begin, cinter.end - bintr.end, NOTES[cinter.data.note]
                ))
        else:
            dotexon.write('{}\t{}\t{}-{}[{}]|{}[{}]\tNovel\t-\t-\t-\t-\n'.format(
                cnt + 1, exon.chrom, cinter.begin, cinter.end - 1, cinter.data.strand, exon.transcript_id,
                cinter.data.transcriptIds[exon.transcript_id]
            ))
    dotexon.close()


def write_itrack(input_file, gtf_chrom_dict, gtf_introns):
    """Writes .itracking files
       num is the number of the current query file"""

    # Create .itracking file
    itrack = open(input_file + '.itracking', 'w')
    if input_file == 'ref':
        interval_best_matches = MatchingDicts.ref_best_matches
        interval_matches = MatchingDicts.ref_matches
        ref_matches = MatchingDicts.interval_matches
        # Write header
        itrack.write("ID\tChromosome\tReference(Coordinates[strand]|Transcript[exon_number])\t"
                    "Query(Ref_Coordinates[strand]|Transcript[exon_number])\tType\t(Start,End)\t[Reference_Structure]"
                    "\t[Query_Structure]\tNotes\n")
    else:
        interval_best_matches = MatchingDicts.interval_best_matches
        interval_matches = MatchingDicts.interval_matches
        ref_matches = MatchingDicts.ref_matches
        # Write header
        itrack.write("ID\tChromosome\tQuery(Coordinates[strand]|Transcript[exon_number])\t"
                    "Reference(Ref_Coordinates[strand]|Transcript[exon_number])\tType\t(Start,End)\t[Query_Structure]"
                    "\t[Reference_Structure]\tNotes\n")
    # Index for each row
    cnt = 0
    # To hold reported intervals not to report them again
    reported = set()
    for chrom in gtf_chrom_dict:
        for strand in ['.', '+', '-']:
            if strand not in gtf_chrom_dict[chrom]:
                continue
            for interval in sorted(gtf_chrom_dict[chrom][strand][0]):
                if interval in reported:
                    continue
                cnt += 1
                if interval in interval_matches:
                    qset = {interval}
                    reported.add(interval)
                    rset = set()
                    # Get the first reference best match's strand
                    gtf_st = six.next(six.iterkeys(interval_best_matches[interval])).data.strand

                    if interval.data.strand != strand:
                        print('SOMETHING IS WRONG STRANDS DOESN\'t MATCH Code: 1')

                    # Get current strand
                    quer_st = interval.data.strand

                    # Add note for the unknown strand if found on both reference strands
                    if len(interval_best_matches[interval]) > 1 and interval.data.note == 7:
                        note = '(A best match on both strands was found only used {} (arbitrary decision))'.format(
                            gtf_st)
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
                                if rinterval.data.strand != gtf_st and rinterval.data.strand != '.':
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
                    for b, e, s in gtf_introns[interval]:
                        if s == strand or strand == '.':
                            lvl = 'Intronic'
                    if interval.data.strand != strand:
                        print('SOMETHING IS WRONG STRANDS DOESN\'t MATCH Code: 2')
                    itrack.write("id_{}\t{}\t{}-{}[{}]|{}\t-\t".format(
                        cnt, chrom, interval.begin, interval.end - 1, interval.data.strand,
                        '|'.join(['{}[{}]'.format(k, v) for k, v in interval.data.transcriptIds.items()])
                    ))
                    itrack.write("SQNR\t-\t-\t-\t[{}]\n".format(lvl))
