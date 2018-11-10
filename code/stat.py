import os
import numpy as np


def get_gap_overlap(type1, type2):
    if type1 == 'Single' or type1 == 'Last':
        first = '3'
    else:
        first = 'intron'

    if type2 == 'Single' or type2 == 'First':
        second = '5'
    else:
        second = 'intron'

    if first == '3' and second == '5':
        return '35prime'
    elif first == second:
        return 'intronic'
    else:
        return first + second


# def get_type(rtype, qtype):
#     if rtype == 'Single':
#         if qtype == 'Single':
#             return ['both', 'both']
#         elif qtype == 'First':
#             return ['both', 'ref']
#         elif qtype == 'Last':
#             return ['ref', 'both']
#         else:
#             return ['ref', 'ref']
#     elif rtype == 'First':
#         if qtype == 'Single':
#             return ['both', 'query']
#         elif qtype == 'First':
#             return ['both', 'intronic']
#         elif qtype == 'Last':
#             return ['ref', 'query']
#         else:
#             return ['ref', 'intronic']
#     elif rtype == 'Last':
#         if qtype == 'Single':
#             return ['query', 'both']
#         elif qtype == 'First':
#             return ['query', 'ref']
#         elif qtype == 'Last':
#             return ['intronic', 'both']
#         else:
#             return ['intronic', 'ref']
#     else:
#         if qtype == 'Single':
#             return ['query', 'query']
#         elif qtype == 'First':
#             return ['query', 'intronic']
#         elif qtype == 'Last':
#             return ['intronic', 'query']
#         else:
#             return ['intronic', 'intronic']

stats = dict()
files = os.listdir(".")
for _file in files:
    allrefints = set()
    allquerints = set()
    if _file.endswith('itracking'):
        print _file
        stats[_file] = {'SQNR': {'count_t': list(), 'count_u': list(), '5_prime_utr': [0, 0, 0, 0], '3_prime_utr': [0, 0, 0, 0],
                                 'ulength': [set(), set(), set()], 'intronic_edge': [0, 0, 0, 0],
                                 'Intronic': {'count_t': list(), 'count_u': list(), 'ulength': [set(), set(), set()]},
                                 'Intergenic': {'count_t': list(), 'count_u': list(), 'ulength': [set(), set(), set()]}
                                 },
                        'SQSR': {'count_t': list(), 'count_u': list(), '5_prime_utr': [0, 0, 0, 0], '3_prime_utr': [0, 0, 0, 0],
                                 'ulength': [set(), set(), set()], 'intronic_edge': [0, 0, 0, 0], '3_prime_ch': [list(), list(), 0],
                                 '5_prime_ch': [list(), list(), 0], 'intronic_ch': [list(), list(), 0]
                                 },
                        'SQMR': {'count_t': list(), 'count_u': list(), '5_prime_utr': [0, 0, 0, 0], '3_prime_utr': [0, 0, 0, 0],
                                 'ulength': [set(), set(), set()], 'intronic_edge': [0, 0, 0, 0], '3_prime_ch': [list(), list(), 0],
                                 '5_prime_ch': [list(), list(), 0], 'intronic_ch': [list(), list(), 0],
                                 'rgapsc': 0, 'rolapsc': 0,
                                 'rgaps': {'intronic': list(), '3intron': list(), 'intron5': list(), '35prime': list()},
                                 'rolaps': {'intronic': list(), '3intron': list(), 'intron5': list(), '35prime': list()}
                                 },
                        'MQSR': {'count_t': list(), 'count_u': list(), '5_prime_utr': [0, 0, 0, 0], '3_prime_utr': [0, 0, 0, 0],
                                 'ulength': [set(), set(), set()], 'intronic_edge': [0, 0, 0, 0], '3_prime_ch': [list(), list(), 0],
                                 '5_prime_ch': [list(), list(), 0], 'intronic_ch': [list(), list(), 0],
                                 'qgapsc': 0, 'qolapsc': 0,
                                 'qgaps': {'intronic': list(), '3intron': list(), 'intron5': list(), '35prime': list()},
                                 'qolaps': {'intronic': list(), '3intron': list(), 'intron5': list(), '35prime': list()}
                                 },
                        'MQMR': {'count_t': list(), 'count_u': list(), '5_prime_utr': [0, 0, 0, 0],
                                 '3_prime_utr': [0, 0, 0, 0], 'ulength': [set(), set(), set()],
                                 'intronic_edge': [0, 0, 0, 0], '3_prime_ch': [list(), list(), 0],
                                 '5_prime_ch': [list(), list(), 0],'intronic_ch': [list(), list(), 0],
                                 'qgapsc': 0, 'qolapsc': 0, 'rgapsc': 0, 'rolapsc': 0,
                                 'rgaps': {'intronic': list(), '3intron': list(), 'intron5': list(), '35prime': list()},
                                 'rolaps': {'intronic': list(), '3intron': list(), 'intron5': list(), '35prime': list()},
                                 'qgaps': {'intronic': list(), '3intron': list(), 'intron5': list(), '35prime': list()},
                                 'qolaps': {'intronic': list(), '3intron': list(), 'intron5': list(), '35prime': list()}
                                 }, 'cnt': 0,
                        'count_t': list(), 'count_u': list(), '5_prime_utr': [0, 0, 0, 0], '3_prime_utr': [0, 0, 0, 0],
                        'ulength': [set(), set(), set()], 'intronic_edge': [0, 0, 0, 0], '3_prime_ch': [list(), list(), 0],
                        '5_prime_ch': [list(), list(), 0], 'intronic_ch': [list(), list(), 0],
                        }
        with open(_file, 'r') as _f:
            for line_index, line in enumerate(_f.readlines()):
                if line_index > 0:
                    columns = line.replace('\n', '').split("\t")
                    _type = columns[4]
                    query_intervals = columns[2].split(', ')
                    reference_intervals = columns[3].split(', ')
                    begin, end = None, None
                    query_stat, ref_stat = None, None
                    if _type != "SQNR":
                        begin, end = columns[5].replace('(', '').replace(')', '').split(',')
                        begin = int(begin)
                        end = int(end)
                        if 'MQ' in _type:
                            query_stat = columns[6].replace('[', '').replace(']', '')
                            query_stat = query_stat.split(' ')
                        if 'MR' in _type:
                            ref_stat = columns[7].replace('[', '').replace(']', '')
                            ref_stat = ref_stat.split(' ')

                    qtypes = list()
                    for single_interval in query_intervals:
                        if single_interval not in allquerints:
                            allquerints.add(single_interval)
                            utr_update = True
                        else:
                            utr_update = False
                            print 'twice', single_interval
                        interval = single_interval.split(']|')
                        trans = interval[1:]
                        interval = interval[0]
                        if _type == 'SQNR':
                            interval = interval.split('][')
                            lvl = interval[1]
                            interval = interval[0]
                        interval, strand = interval[0:-1].replace('[', ''), interval[-1]
                        query_transcripts = dict()
                        qitype = set()
                        for i, item in enumerate(trans):
                            k, v = item.split('[')
                            v = v.replace(']', '')
                            if k in query_transcripts:
                                if query_transcripts[k] == v:
                                    print "Error dual transcript", k
                                    exit()
                            query_transcripts[k] = v
                            qitype.add(v)

                            if utr_update:
                                if v == 'Single':
                                    stats[_file]['5_prime_utr'][0] += 1
                                    stats[_file]['3_prime_utr'][0] += 1
                                    stats[_file][_type]['5_prime_utr'][0] += 1
                                    stats[_file][_type]['3_prime_utr'][0] += 1

                                elif v == 'First':
                                    stats[_file]['5_prime_utr'][0] += 1
                                    stats[_file]['intronic_edge'][0] += 1
                                    stats[_file][_type]['5_prime_utr'][0] += 1
                                    stats[_file][_type]['intronic_edge'][0] += 1
                                elif v == 'Last':
                                    stats[_file]['3_prime_utr'][0] += 1
                                    stats[_file]['intronic_edge'][0] += 1
                                    stats[_file][_type]['3_prime_utr'][0] += 1
                                    stats[_file][_type]['intronic_edge'][0] += 1
                                else:
                                    stats[_file]['intronic_edge'][0] += 2
                                    stats[_file][_type]['intronic_edge'][0] += 2

                        if 'Single' in qitype:
                            qtypes.append('Single')
                            if utr_update:
                                stats[_file]['5_prime_utr'][1] += 1
                                stats[_file]['3_prime_utr'][1] += 1
                                stats[_file][_type]['5_prime_utr'][1] += 1
                                stats[_file][_type]['3_prime_utr'][1] += 1
                        elif 'First' in qitype:
                            qtypes.append('First')
                            if utr_update:
                                stats[_file]['5_prime_utr'][1] += 1
                                stats[_file]['intronic_edge'][1] += 1
                                stats[_file][_type]['5_prime_utr'][1] += 1
                                stats[_file][_type]['intronic_edge'][1] += 1
                        elif 'Last' in qitype:
                            qtypes.append('Last')
                            if utr_update:
                                stats[_file]['3_prime_utr'][1] += 1
                                stats[_file]['intronic_edge'][1] += 1
                                stats[_file][_type]['3_prime_utr'][1] += 1
                                stats[_file][_type]['intronic_edge'][1] += 1
                        else:
                            qtypes.append('Mid')
                            if utr_update:
                                stats[_file]['intronic_edge'][1] += 2
                                stats[_file][_type]['intronic_edge'][1] += 2

                        if len(qtypes) > 1:
                            go_type = get_gap_overlap(qtypes[-2], qtypes[-1])
                            # print query_stat, _type, '\n', line_index
                            qstat = int(query_stat[len(qtypes) - 2])
                            if qstat < 0:
                                stats[_file][_type]['qgaps'][go_type].append(qstat)
                                stats[_file][_type]['qgapsc'] += 1
                            else:
                                stats[_file][_type]['qolaps'][go_type].append(qstat)
                                stats[_file][_type]['qolapsc'] += 1

                        interval = interval.split('-')
                        size = int(interval[1]) - int(interval[0]) + 1

                        if len(query_transcripts) == 0:
                            print line_index, _type
                            exit()
                        stats[_file][_type]['count_t'] += [size] * len(query_transcripts)
                        stats[_file]['count_t'] += [size] * len(query_transcripts)

                        stats[_file][_type]['count_u'].append(size)
                        stats[_file]['count_u'].append(size)
                        if _type == 'SQNR':
                            if not lvl:
                                print 'SQNR without level at', interval
                                exit()
                            else:
                                stats[_file][_type][lvl]['count_t'] += [size] * len(query_transcripts)
                                stats[_file][_type][lvl]['count_u'].append(size)
                        intersect = list()
                        if strand == '+':
                            usl = 0
                        elif strand == '-':
                            usl = 1
                        else:
                            usl = 2
                        for unique in range(int(interval[0]), int(interval[1]) + 1):
                            if unique not in stats[_file]['ulength'][usl]:
                                stats[_file]['ulength'][usl].add(unique)
                                stats[_file][_type]['ulength'][usl].add(unique)
                                if _type == 'SQNR':
                                    if not lvl:
                                        print 'SQNR without level at', interval
                                        exit()
                                    else:
                                        stats[_file][_type][lvl]['ulength'][usl].add(unique)

                    rtypes = list()
                    mxr = 0
                    for single_interval in reference_intervals:
                        if single_interval not in allrefints:
                            allrefints.add(single_interval)
                            utr_update = True
                        else:
                            utr_update = False
                        if single_interval != '-':
                            interval = single_interval.split(']|')
                            ref_transcripts = dict()
                            ritype = set()
                            for i, item in enumerate(interval):
                                if i > 0:
                                    item = item.split('[')
                                    k, v = item
                                    v = v.replace(']', '')
                                    if k in ref_transcripts:
                                        if ref_transcripts[k] == v:
                                            print "Error dual transcript", k
                                            exit()
                                    ref_transcripts[k] = v
                                    ritype.add(v)

                                    if utr_update:
                                        if v == 'Single':
                                            stats[_file]['5_prime_utr'][2] += 1
                                            stats[_file]['3_prime_utr'][2] += 1
                                            stats[_file][_type]['5_prime_utr'][2] += 1
                                            stats[_file][_type]['3_prime_utr'][2] += 1

                                        elif v == 'First':
                                            stats[_file]['5_prime_utr'][2] += 1
                                            stats[_file]['intronic_edge'][2] += 1
                                            stats[_file][_type]['5_prime_utr'][2] += 1
                                            stats[_file][_type]['intronic_edge'][2] += 1
                                        elif v == 'Last':
                                            stats[_file]['3_prime_utr'][2] += 1
                                            stats[_file]['intronic_edge'][2] += 1
                                            stats[_file][_type]['3_prime_utr'][2] += 1
                                            stats[_file][_type]['intronic_edge'][2] += 1
                                        else:
                                            stats[_file]['intronic_edge'][2] += 2
                                            stats[_file][_type]['intronic_edge'][2] += 2

                            if 'Single' in ritype:
                                rtypes.append('Single')
                                if utr_update:
                                    stats[_file]['5_prime_utr'][3] += 1
                                    stats[_file]['3_prime_utr'][3] += 1
                                    stats[_file][_type]['5_prime_utr'][3] += 1
                                    stats[_file][_type]['3_prime_utr'][3] += 1
                            elif 'First' in ritype:
                                rtypes.append('First')
                                if utr_update:
                                    stats[_file]['5_prime_utr'][3] += 1
                                    stats[_file]['intronic_edge'][3] += 1
                                    stats[_file][_type]['5_prime_utr'][3] += 1
                                    stats[_file][_type]['intronic_edge'][3] += 1
                            elif 'Last' in ritype:
                                rtypes.append('Last')
                                if utr_update:
                                    stats[_file]['3_prime_utr'][3] += 1
                                    stats[_file]['intronic_edge'][3] += 1
                                    stats[_file][_type]['3_prime_utr'][3] += 1
                                    stats[_file][_type]['intronic_edge'][3] += 1
                            else:
                                rtypes.append('Mid')
                                if utr_update:
                                    stats[_file]['intronic_edge'][3] += 2
                                    stats[_file][_type]['intronic_edge'][3] += 2

                            if len(rtypes) > 1:
                                go_type = get_gap_overlap(rtypes[-2], rtypes[-1])
                                rstat = int(ref_stat[len(rtypes) - 2])
                                if rstat < 0:
                                    stats[_file][_type]['rgaps'][go_type].append(rstat)
                                    stats[_file][_type]['rgapsc'] += 1
                                else:
                                    stats[_file][_type]['rolaps'][go_type].append(rstat)
                                    stats[_file][_type]['rolapsc'] += 1

                            interval = interval[0][2:-2]
                            interval = interval.split('-')
                            if int(interval[1]) > mxr:
                                last_ref = rtypes[-1]
                                mxr = int(interval[1])

                    if _type != 'SQNR':

                        if begin > 0:
                            x = 0
                        elif begin < 0:
                            x = 1
                        else:
                            x = 2

                        if rtypes[0] == 'Single' or rtypes[0] == 'First':
                            change = '5_prime_ch'
                        else:
                            change = 'intronic_ch'

                        if x != 2:
                            stats[_file][change][x].append(abs(begin))
                            stats[_file][_type][change][x].append(abs(begin))
                        else:
                            stats[_file][change][x] += 1
                            stats[_file][_type][change][x] += 1

                        if end > 0:
                            x = 0
                        elif end < 0:
                            x = 1
                        else:
                            x = 2

                        if last_ref == 'Single' or last_ref == 'Last':
                            change = '3_prime_ch'
                        else:
                            change = 'intronic_ch'

                        if x != 2:
                            stats[_file][change][x].append(abs(end))
                            stats[_file][_type][change][x].append(abs(end))
                        else:
                            stats[_file][change][x] += 1
                            stats[_file][_type][change][x] += 1
        for uil in range(3):
            stats[_file]['ulength'][uil] = len(stats[_file]['ulength'][uil])
            for kk in ['SQNR', 'SQSR', 'SQMR', 'MQSR', 'MQMR']:
                stats[_file][kk]['ulength'][uil] = len(stats[_file][kk]['ulength'][uil])
                if kk == 'SQNR':
                    stats[_file][kk]['Intronic']['ulength'][uil] = len(stats[_file][kk]['Intronic']['ulength'][uil])
                    stats[_file][kk]['Intergenic']['ulength'][uil] = len(stats[_file][kk]['Intergenic']['ulength'][uil])

print stats['CQUERY_SAMPLE_TEST.gtf.db.itracking']['5_prime_utr'][0], stats['CQUERY_SAMPLE_TEST.gtf.db.itracking']['3_prime_utr'][0]
print stats['QUERY_SAMPLE_TEST.gtf.db.itracking']['5_prime_utr'][0], stats['QUERY_SAMPLE_TEST.gtf.db.itracking']['3_prime_utr'][0]
with open('stats.txt', 'w') as write:
    for asm in ['stringtie', 'scallop', 'trinity']:
        write.write('\t' * 4 + asm + '\n')
        write.write('\t' * 2 + 'poly' + '\t' * 3 + 'ribo' + '\n')
        write.write('\t' * 1 + 'brain' + '\t' * 1 + 'tumor')
        write.write('\t' * 1 + 'brain' + '\t' * 1 + 'tumor' + '\n')
        lines = [''] * 150
        lines[0] = 'Total_Exons(count|total|mean|std)\t'
        lines[1] = 'Unique_Exons(count|total|mean|std)\t'
        lines[2] = 'Transcripts(count)\t'
        lines[3] = 'Intronic_Edges(count)\t'

        xx = 4
        for stat in ['SQNR', 'SQSR', 'SQMR', 'MQSR', 'MQMR']:
            for tu in ['Total', 'Unique']:
                lines[xx] = '{}_{}_Exons(count(%)|total|mean|std)\t'.format(stat, tu)
                xx += 1

                if stat == 'SQNR':
                    for sstat in ['Intronic', 'Intergenic']:
                        lines[xx] = '{}_{}_{}_Exons(count(%)|total|mean|std)\t'.format(stat, sstat, tu)
                        xx += 1
                elif tu == 'Total':
                    for mutr in ['5_prime', '3_prime', 'Intronic_Egdes']:
                        if mutr != 'Intronic_Edges':
                            utr = '_UTR'
                        else:
                            utr = ''
                        lines[xx] = '{}_{}_Matching_{}{}(count(%))\t'.format(stat, tu, mutr, utr)
                        xx += 1

            if stat != 'SQNR':
                for mutr in ['5_prime', '3_prime', 'Intronic_Egdes']:
                    if mutr != 'Intronic_Edges':
                        utr = '_UTR'
                    else:
                        utr = ''
                    lines[xx] = '{}_{}_Matching_{}{}(count(%))\t'.format(stat, tu, mutr, utr)
                    xx += 1
                    for k in ['Positive', 'Negative', 'No_Change']:
                        lines[xx] = '{}_{}_{}(count(%)|total|mean|std)\t'.format(stat, mutr, k)
                        xx += 1

            if 'MQ' in stat:
                for go in ['qgap', 'qoverlap']:
                    lines[xx] = '{}_{}s\t'.format(stat, go)
                    xx += 1
                    for sstat in ['intronic', '3intron', 'intron5', '35prime']:
                        lines[xx] = '{}_{}_{}(count(%)|total|mean|std)\t'.format(stat, go, sstat)
                        xx += 1

            if 'MR' in stat:
                for go in ['rgap', 'roverlap']:
                    lines[xx] = '{}_{}s\t'.format(stat, go)
                    xx += 1
                    for sstat in ['intronic', '3intron', 'intron5', '35prime']:
                        lines[xx] = '{}_{}_{}(count(%)|total|mean|std)\t'.format(stat, go, sstat)
                        xx += 1
        print xx, 'sad'

        for lib in ['poly', 'ribo']:
            for sample in ['brain', 'tumor']:
                for _file in stats:
                    if asm in _file and lib in _file and sample in _file:

                        val = stats[_file]['count_t']
                        val = np.array(val)
                        lines[0] += '{}|{}|{}|{}\t'.format(len(val), val.sum(),
                                                           round(val.mean(), 1), round(val.std(), 1))
                        val = stats[_file]['count_u']
                        val = np.array(val)
                        val2 = sum(stats[_file]['ulength'])
                        lines[1] += '{}|{}|{}|{}\t'.format(len(val), val2,
                                                           round(val.mean(), 1), round(val.std(), 1))
                        val, val2 = stats[_file]['5_prime_utr'][0], stats[_file]['3_prime_utr'][0]
                        if val != val2:
                            print 'UTRs don\'t match', val, val2
                        lines[2] += '{}\t'.format(val)
                        val = stats[_file]['intronic_edge'][0]
                        lines[3] += '{}\t'.format(val)

                        xx = 4
                        for stat in ['SQNR', 'SQSR', 'SQMR', 'MQSR', 'MQMR']:
                            for tt in ['count_t', 'count_u']:
                                val = stats[_file][stat][tt]
                                val = np.array(val)
                                if tt == 'count_u':
                                    val2 = sum(stats[_file][stat]['ulength'])
                                else:
                                    val2 = val.sum()
                                per = len(val) * 100.0 / len(stats[_file][tt])
                                lines[xx] += '{}({})|{}|{}|{}\t'.format(len(val), round(per, 1), val2,
                                                                        round(val.mean(), 1), round(val.std(), 1))
                                xx += 1

                                if stat == 'SQNR':
                                    for sstat in ['Intronic', 'Intergenic']:
                                        val = stats[_file][stat][sstat][tt]
                                        val = np.array(val)
                                        per = len(val) * 100.0 / len(stats[_file][stat][tt])
                                        lines[xx] += '{}({})|{}|{}|{}\t'.format(len(val), round(per, 1), val.sum(),
                                                                            round(val.mean(), 1), round(val.std(), 1))
                                        xx += 1
                                elif tt == 'count_t':
                                    for mutr in ['5_prime_utr', '3_prime_utr', 'intronic_edge']:
                                        val = stats[_file][stat][mutr][2]
                                        per = val * 100.0 / stats[_file][mutr][2]
                                        lines[xx] += '{}({})\t'.format(val, round(per, 1))
                                        xx += 1

                            if stat != 'SQNR':
                                for mutr in ['5_prime_utr', '3_prime_utr', 'intronic_edge']:
                                    val2 = stats[_file][stat][mutr][3]
                                    per = val2 * 100.0 / stats[_file][mutr][3]
                                    lines[xx] += '{}({})\t'.format(val2, round(per, 1))
                                    xx += 1
                                    mutr = mutr.replace('utr', 'ch').replace('edge', 'ch')
                                    for k in range(2):
                                        val = stats[_file][stat][mutr][k]
                                        val = np.array(val)
                                        val3 = len(val) * 100.0 / val2
                                        lines[xx] += '{}({})|{}|{}|{}\t'.format(len(val), round(val3, 1), val.sum(),
                                                                                round(val.mean(), 1),
                                                                                round(val.std(), 1))
                                        xx += 1
                                    val = stats[_file][stat][mutr][2]
                                    val3 = val * 100.0 / val2
                                    lines[xx] += '{}({})\t'.format(val, round(val3, 1))
                                    xx += 1

                            if 'MQ' in stat:
                                for go in ['qgaps', 'qolaps']:
                                    val2 = stats[_file][stat][go + 'c']
                                    lines[xx] += '{}\t'.format(val2)
                                    xx += 1
                                    for sstat in ['intronic', '3intron', 'intron5', '35prime']:
                                        val = stats[_file][stat][go][sstat]
                                        val = np.array(val)
                                        per = len(val) * 100.0 / val2
                                        lines[xx] += '{}({})|{}|{}|{}\t'.format(len(val), round(per, 1), val.sum(),
                                                                            round(val.mean(), 1), round(val.std(), 1))
                                        xx += 1
                            if 'MR' in stat:
                                for go in ['rgaps', 'rolaps']:
                                    val2 = stats[_file][stat][go + 'c']
                                    lines[xx] += '{}\t'.format(val2)
                                    xx += 1
                                    for sstat in ['intronic', '3intron', 'intron5', '35prime']:
                                        val = stats[_file][stat][go][sstat]
                                        val = np.array(val)
                                        per = len(val) * 100.0 / val2
                                        lines[xx] += '{}({})|{}|{}|{}\t'.format(len(val), round(per, 1), val.sum(),
                                                                            round(val.mean(), 1), round(val.std(), 1))
                                        xx += 1
        print xx, 'here'

        for line in lines:
            write.write(line + '\n')

        write.write('\n')
        write.write('\n')
