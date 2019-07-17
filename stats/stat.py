import os
import numpy as np
import time
from intervaltree import IntervalTree, Interval
import openpyxl

np.seterr(divide='ignore', invalid='ignore')

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


stats = {}
rnd = 1
files = os.listdir(".")
for _file in files:
    allrefints = set()
    allquerints = set()
    time1 = time.time()
    if _file.endswith('itracking'):

        stats[_file] = {'SQNR': {'count_t': [], 'count_u': [], '5_prime_utr': [0, 0, 0, []],
                                 '3_prime_utr': [0, 0, 0, []], 'intronic_edge': [0, 0, 0, []],
                                 '5_prime_ch_u': [0, 0, 0], '3_prime_ch_u': [0, 0, 0], 'intronic_ch_u': [0, 0, 0],
                                 '5_prime_utr_u': 0, '3_prime_utr_u': 0, 'intronic_edge_u': 0,
                                 'ulength': {'+': IntervalTree(), '-': IntervalTree(), '.': IntervalTree()},
                                 'Intronic': {
                                     'count_t': [], 'count_u': [],
                                     'ulength': {'+': IntervalTree(), '-': IntervalTree(), '.': IntervalTree()}},
                                 'Intergenic': {
                                     'count_t': [], 'count_u': [],
                                     'ulength': {'+': IntervalTree(), '-': IntervalTree(), '.': IntervalTree()}}
                                 },
                        'SQSR': {'count_t': [], '5_prime_utr': [0, 0, 0, []], '3_prime_utr': [0, 0, 0, []],
                                 'count_u': [], 'intronic_edge': [0, 0, 0, []], '3_prime_ch': [[], [], []],
                                 'ulength': {'+': IntervalTree(), '-': IntervalTree(), '.': IntervalTree()},
                                 '5_prime_ch': [[], [], []], 'intronic_ch': [[], [], []],
                                 '5_prime_ch_u': [0, 0, 0], '3_prime_ch_u': [0, 0, 0], 'intronic_ch_u': [0, 0, 0],
                                 '5_prime_utr_u': 0, '3_prime_utr_u': 0, 'intronic_edge_u': 0
                                 },
                        'SQMR': {'count_t': [], '5_prime_utr': [0, 0, 0, []], '3_prime_utr': [0, 0, 0, []],
                                 'ulength': {'+': IntervalTree(), '-': IntervalTree(), '.': IntervalTree()},
                                 'intronic_edge': [0, 0, 0, []], '3_prime_ch': [[], [], []], '5_prime_ch': [[], [], []],
                                 '5_prime_ch_u': [0, 0, 0], '3_prime_ch_u': [0, 0, 0], 'intronic_ch_u': [0, 0, 0],
                                 '5_prime_utr_u': 0, '3_prime_utr_u': 0, 'intronic_edge_u': 0,
                                 'intronic_ch': [[], [], []], 'rgapsc': [], 'rolapsc': [], 'count_u': [],
                                 'rgaps': {'intronic': [], '3intron': [], 'intron5': [], '35prime': []},
                                 'rolaps': {'intronic': [], '3intron': [], 'intron5': [], '35prime': []}
                                 },
                        'MQSR': {'count_t': [], '5_prime_utr': [0, 0, 0, []], '3_prime_utr': [0, 0, 0, []],
                                 'ulength': {'+': IntervalTree(), '-': IntervalTree(), '.': IntervalTree()},
                                 'intronic_edge': [0, 0, 0, []], '3_prime_ch': [[], [], []], '5_prime_ch': [[], [], []],
                                 '5_prime_ch_u': [0, 0, 0], '3_prime_ch_u': [0, 0, 0], 'intronic_ch_u': [0, 0, 0],
                                 '5_prime_utr_u': 0, '3_prime_utr_u': 0, 'intronic_edge_u': 0,
                                 'intronic_ch': [[], [], []], 'qgapsc': [], 'qolapsc': [], 'count_u': [],
                                 'qgaps': {'intronic': [], '3intron': [], 'intron5': [], '35prime': []},
                                 'qolaps': {'intronic': [], '3intron': [], 'intron5': [], '35prime': []}
                                 },
                        'MQMR': {'count_t': [], '5_prime_utr': [0, 0, 0, []], '3_prime_utr': [0, 0, 0, []],
                                 'ulength': {'+': IntervalTree(), '-': IntervalTree(), '.': IntervalTree()},
                                 'intronic_edge': [0, 0, 0, []], '3_prime_ch': [[], [], []], '5_prime_ch': [[], [], []],
                                 '5_prime_ch_u': [0, 0, 0], '3_prime_ch_u': [0, 0, 0], 'intronic_ch_u': [0, 0, 0],
                                 '5_prime_utr_u': 0, '3_prime_utr_u': 0, 'intronic_edge_u': 0,
                                 'intronic_ch': [[], [], []], 'qgapsc': [], 'qolapsc': [], 'rgapsc': [], 'rolapsc': [],
                                 'rgaps': {'intronic': [], '3intron': [], 'intron5': [], '35prime': []}, 'count_u': [],
                                 'rolaps': {'intronic': [], '3intron': [], 'intron5': [], '35prime': []},
                                 'qgaps': {'intronic': [], '3intron': [], 'intron5': [], '35prime': []},
                                 'qolaps': {'intronic': [], '3intron': [], 'intron5': [], '35prime': []}
                                 }, 'cnt': 0, 'count_t': [], '5_prime_utr': [0, 0, 0, []], '3_prime_utr': [0, 0, 0, []],
                        'count_u': [], 'ulength': {'+': IntervalTree(), '-': IntervalTree(), '.': IntervalTree()},
                        'intronic_edge': [0, 0, 0, []], '3_prime_ch': [[], [], []], '5_prime_ch': [[], [], []],
                        'intronic_ch': [[], [], []], '5_prime_ch_u': [0, 0, 0], '3_prime_ch_u': [0, 0, 0],
                        'intronic_ch_u': [0, 0, 0], '5_prime_utr_u': 0, '3_prime_utr_u': 0, 'intronic_edge_u': 0
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
                    qtypes = []
                    for single_interval in query_intervals:
                        if single_interval not in allquerints:
                            allquerints.add(single_interval)
                            utr_update = True
                        else:
                            utr_update = False
                            print('twice', single_interval)
                        interval = single_interval.split(']|')
                        trans = interval[1:]
                        interval = interval[0]
                        if _type == 'SQNR':
                            # interval = interval.split('][')
                            # lvl = interval[1]
                            interval = interval.split(']')
                            lvl = columns[8].replace('[', '').replace(']', '')
                            interval = interval[0]
                        interval, strand = interval[0:-1].replace('[', ''), interval[-1]
                        query_transcripts = {}
                        qitype = set()
                        for i, item in enumerate(trans):
                            k, v = item.split('[')
                            v = v.replace(']', '')
                            if k in query_transcripts:
                                if query_transcripts[k] == v:
                                    print("Error dual transcript", k)
                                    exit()
                            query_transcripts[k] = v
                            qitype.add(v)

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
                            stats[_file]['5_prime_utr'][1] += 1
                            stats[_file]['3_prime_utr'][1] += 1
                            stats[_file][_type]['5_prime_utr'][1] += 1
                            stats[_file][_type]['3_prime_utr'][1] += 1

                        elif 'First' in qitype:
                            qtypes.append('First')
                            stats[_file]['5_prime_utr'][1] += 1
                            stats[_file]['intronic_edge'][1] += 1
                            stats[_file][_type]['5_prime_utr'][1] += 1
                            stats[_file][_type]['intronic_edge'][1] += 1

                        elif 'Last' in qitype:
                            qtypes.append('Last')
                            stats[_file]['3_prime_utr'][1] += 1
                            stats[_file]['intronic_edge'][1] += 1
                            stats[_file][_type]['3_prime_utr'][1] += 1
                            stats[_file][_type]['intronic_edge'][1] += 1
                        else:
                            qtypes.append('Mid')
                            stats[_file]['intronic_edge'][1] += 2
                            stats[_file][_type]['intronic_edge'][1] += 2

                        if len(qtypes) > 1:
                            go_type = get_gap_overlap(qtypes[-2], qtypes[-1])
                            # print(query_stat, _type, '\n', line_index)
                            qstat = int(query_stat[len(qtypes) - 2])
                            if qstat < 0:
                                stats[_file][_type]['qgaps'][go_type].append(abs(qstat))
                                stats[_file][_type]['qgapsc'].append(abs(qstat))
                            else:
                                stats[_file][_type]['qolaps'][go_type].append(qstat)
                                stats[_file][_type]['qolapsc'].append(qstat)

                        interval = interval.split('-')
                        interval[1], interval[0] = int(interval[1]), int(interval[0])
                        size = interval[1] - interval[0] + 1

                        if len(query_transcripts) == 0:
                            print(line_index, _type)
                            exit()
                        stats[_file][_type]['count_t'] += [size] * len(query_transcripts)
                        stats[_file]['count_t'] += [size] * len(query_transcripts)

                        if _type == 'SQNR':
                            if not lvl:
                                print('SQNR without level at', interval)
                                exit()
                            else:
                                stats[_file][_type][lvl]['count_t'] += [size] * len(query_transcripts)

                        mn = interval[0]
                        mx = interval[1] + 1
                        for b, e, d in stats[_file]['ulength'][strand][mn - 1:mx]:
                            if b < mn:
                                mn = b
                            if e > mx:
                                mx = e
                            stats[_file]['ulength'][strand].removei(b, e, d)
                            stats[_file][_type]['ulength'][strand].discardi(b, e, d)
                            try:
                                stats[_file]['count_u'].remove(e - b)
                                stats[_file][_type]['count_u'].remove(e - b)
                            except ValueError:
                                pass
                            if _type == 'SQNR':
                                if not lvl:
                                    print('SQNR without level at', interval)
                                    exit()
                                else:
                                    stats[_file][_type][lvl]['ulength'][strand].discardi(b, e, d)
                                    try:
                                        stats[_file][_type][lvl]['count_u'].remove(e - b)
                                    except ValueError:
                                        pass
                        stats[_file]['ulength'][strand].addi(mn, mx, i)
                        stats[_file][_type]['ulength'][strand].addi(mn, mx, i)
                        stats[_file]['count_u'].append(mx - mn)
                        stats[_file][_type]['count_u'].append(mx - mn)
                        if _type == 'SQNR':
                            stats[_file][_type][lvl]['ulength'][strand].addi(mn, mx, i)
                            stats[_file][_type][lvl]['count_u'].append(mx - mn)
                    rtypes = []
                    utr_updates = []
                    mxr = 0
                    for single_interval in reference_intervals:
                        if single_interval != '-':
                            if single_interval not in allrefints:
                                allrefints.add(single_interval)
                                utr_update = True
                                utr_updates.append(True)
                            else:
                                utr_update = False
                                utr_updates.append(False)
                                # print('twice', single_interval)
                            interval = single_interval.split(']|')
                            ref_transcripts = {}
                            ritype = set()
                            for i, item in enumerate(interval):
                                if i > 0:
                                    item = item.split('[')
                                    k, v = item
                                    v = v.replace(']', '')
                                    if k in ref_transcripts:
                                        if ref_transcripts[k] == v:
                                            print("Error dual transcript", k)
                                            exit()
                                    ref_transcripts[k] = v
                                    ritype.add(v)

                                    if not utr_update:
                                        continue

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
                                    stats[_file]['5_prime_utr_u'] += 1
                                    stats[_file]['3_prime_utr_u'] += 1
                                    stats[_file][_type]['5_prime_utr_u'] += 1
                                    stats[_file][_type]['3_prime_utr_u'] += 1
                                stats[_file]['5_prime_utr'][3].append(abs(begin))
                                stats[_file]['3_prime_utr'][3].append(abs(end))
                                stats[_file][_type]['5_prime_utr'][3].append(abs(begin))
                                stats[_file][_type]['3_prime_utr'][3].append(abs(end))
                            elif 'First' in ritype:
                                rtypes.append('First')
                                if utr_update:
                                    stats[_file]['5_prime_utr_u'] += 1
                                    stats[_file]['intronic_edge_u'] += 1
                                    stats[_file][_type]['5_prime_utr_u'] += 1
                                    stats[_file][_type]['intronic_edge_u'] += 1
                                stats[_file]['5_prime_utr'][3].append(abs(begin))
                                stats[_file]['intronic_edge'][3].append(abs(end))
                                stats[_file][_type]['5_prime_utr'][3].append(abs(begin))
                                stats[_file][_type]['intronic_edge'][3].append(abs(end))

                            elif 'Last' in ritype:
                                rtypes.append('Last')
                                if utr_update:
                                    stats[_file]['3_prime_utr_u'] += 1
                                    stats[_file]['intronic_edge_u'] += 1
                                    stats[_file][_type]['3_prime_utr_u'] += 1
                                    stats[_file][_type]['intronic_edge_u'] += 1
                                stats[_file]['3_prime_utr'][3].append(abs(end))
                                stats[_file]['intronic_edge'][3].append(abs(begin))
                                stats[_file][_type]['3_prime_utr'][3].append(abs(end))
                                stats[_file][_type]['intronic_edge'][3].append(abs(begin))
                            else:
                                rtypes.append('Mid')
                                if utr_update:
                                    stats[_file]['intronic_edge_u'] += 2
                                    stats[_file][_type]['intronic_edge_u'] += 2
                                stats[_file]['intronic_edge'][3] += [abs(begin), abs(end)]
                                stats[_file][_type]['intronic_edge'][3] += [abs(begin), abs(end)]

                            if len(rtypes) > 1:
                                go_type = get_gap_overlap(rtypes[-2], rtypes[-1])
                                rstat = int(ref_stat[len(rtypes) - 2])
                                if rstat < 0:
                                    stats[_file][_type]['rgaps'][go_type].append(abs(rstat))
                                    stats[_file][_type]['rgapsc'].append(abs(rstat))
                                else:
                                    stats[_file][_type]['rolaps'][go_type].append(rstat)
                                    stats[_file][_type]['rolapsc'].append(rstat)

                            interval = interval[0][2:-2]
                            interval = interval.split('-')
                            if int(interval[1]) > mxr:
                                last_ref = rtypes[-1]
                                last_utr_update = utr_updates[-1]
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

                        if utr_updates[0]:
                            stats[_file][change + '_u'][x] += 1
                            stats[_file][_type][change + '_u'][x] += 1

                        stats[_file][change][x].append(abs(begin))
                        stats[_file][_type][change][x].append(abs(begin))

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

                        if last_utr_update:
                            stats[_file][change + '_u'][x] += 1
                            stats[_file][_type][change + '_u'][x] += 1

                        stats[_file][change][x].append(abs(end))
                        stats[_file][_type][change][x].append(abs(end))

        # print(_file, time.time() - time1)

with open('stats.txt', 'w') as write:
    # write = csv.writer(write, delimiter='\t', quotechar=',')
    with open('refstat.txt', 'r') as rr:
        lines = rr.readlines()
        for line in lines:
            write.write('REF_' + line)
    lines = [''] * 3

    # for asm in ['stringtie', 'scallop', 'trinity']:
    #     lines[0] += '\t' * 4 + asm + '\t#'
    #     lines[1] += '\t' * 2 + 'poly' + '\t' * 3 + 'ribo' + '\t#'
    #     lines[2] += '\t' * 1 + 'brain' + '\t' * 1 + 'tumor'
    #     lines[2] += '\t' * 1 + 'brain' + '\t' * 1 + 'tumor' + '\t#'
    # for l in lines:
    #     write.write(l)
    #     write.write('\n')
    lines = [''] * 119
    lines[0] = 'Total_Exons(count|total|mean|std)\t'
    lines[1] = 'Unique_Exons(count|total|mean|std)\t'
    lines[2] = 'Transcripts(count)\t'
    lines[3] = 'Intronic_Edges(count)\t'

    xx = 4
    for stat in ['SQNR', 'SQSR', 'SQMR', 'MQSR', 'MQMR']:
        for tu in ['Total', 'Unique']:
            lines[xx] = '{}_{}_Exons(count(%)|total(%)|mean|std)\t'.format(stat, tu)
            xx += 1

            if stat == 'SQNR':
                for sstat in ['Intronic', 'Intergenic']:
                    lines[xx] = '{}_{}_{}_Exons(count(%)|total(%)|mean|std)\t'.format(stat, sstat, tu)
                    xx += 1
            elif tu == 'Total':
                for mutr in ['5_prime', '3_prime', 'Intronic_Edges']:
                    if mutr != 'Intronic_Edges':
                        utr = '_UTR'
                    else:
                        utr = ''
                    lines[xx] = '{}_{}_Matching_{}{}(count(%)|total(%))\t'.format(stat, tu, mutr, utr)
                    xx += 1

        if stat != 'SQNR':
            for mutr in ['5_prime', '3_prime', 'Intronic_Edges']:
                if mutr != 'Intronic_Edges':
                    utr = '_UTR'
                else:
                    utr = ''
                lines[xx] = '{}_{}_Matching_{}{}(count(%)|total(%))\t'.format(stat, tu, mutr, utr)
                xx += 1
                for k in ['Positive', 'Negative']:
                    lines[xx] = '{}_{}_{}(count(%)|total(%)|mean|std)\t'.format(stat, mutr, k)
                    xx += 1
                lines[xx] = '{}_{}_{}(count(%))\t'.format(stat, mutr, 'No_Change')
                xx += 1

        if 'MQ' in stat:
            for go in ['qgap', 'qoverlap']:
                lines[xx] = '{}_{}s(count|total)\t'.format(stat, go)
                xx += 1
                for sstat in ['intronic', '3intron', 'intron5', '35prime']:
                    lines[xx] = '{}_{}_{}(count(%)|total(%)|mean|std)\t'.format(stat, go, sstat)
                    xx += 1

        if 'MR' in stat:
            for go in ['rgap', 'roverlap']:
                lines[xx] = '{}_{}s(count|total)\t'.format(stat, go)
                xx += 1
                for sstat in ['intronic', '3intron', 'intron5', '35prime']:
                    lines[xx] = '{}_{}_{}(count(%)|total(%)|mean|std)\t'.format(stat, go, sstat)
                    xx += 1

    # for asm in ['stringtie', 'scallop', 'trinity']:
    #     for lib in ['poly', 'ribo']:
    #         for sample in ['brain', 'tumor']:
    for _file in stats:
                    # if asm in _file and lib in _file and sample in _file:

        val = stats[_file]['count_t']
        val = np.array(val)
        lines[0] += '{}|{}|{}|{}\t'.format(len(val), val.sum(),
                                            round(val.mean(), 1), round(val.std(), 1))
        val = stats[_file]['count_u']
        val = np.array(val)
        # val2 = sum(stats[_file]['ulength'])
        lines[1] += '{}|{}|{}|{}\t'.format(len(val), val.sum(),
                                            round(val.mean(), 1), round(val.std(), 1))
        val, val2 = stats[_file]['5_prime_utr'][0], stats[_file]['3_prime_utr'][0]
        if val != val2:
            print('UTRs don\'t match', val, val2)
        lines[2] += '{}\t'.format(val)
        val = stats[_file]['intronic_edge'][0]
        lines[3] += '{}\t'.format(val)

        xx = 4
        for stat in ['SQNR', 'SQSR', 'SQMR', 'MQSR', 'MQMR']:
            for tt in ['count_t', 'count_u']:
                val = stats[_file][stat][tt]
                val = np.array(val)
                # if tt == 'count_u':
                #     val2 = sum(stats[_file][stat]['ulength'])
                # else:
                #     val2 = val.sum()
                perc = len(val) * 100.0 / len(stats[_file][tt])
                pert = sum(val) * 100.0 / sum(stats[_file][tt])
                lines[xx] += '{}({})|{}({})|{}|{}\t'.format(len(val), round(perc, rnd), val.sum(),
                                                            round(pert, rnd), round(val.mean(), rnd),
                                                            round(val.std(), rnd))
                xx += 1

                if stat == 'SQNR':
                    for sstat in ['Intronic', 'Intergenic']:
                        val = stats[_file][stat][sstat][tt]
                        val = np.array(val)
                        perc = len(val) * 100.0 / len(stats[_file][stat][tt])
                        pert = sum(val) * 100.0 / sum(stats[_file][stat][tt])
                        lines[xx] += '{}({})|{}({})|{}|{}\t'.format(len(val), round(perc, rnd),
                                                                    val.sum(), round(pert, rnd),
                                                                    round(val.mean(), rnd),
                                                                    round(val.std(), rnd))
                        xx += 1
                elif tt == 'count_t':
                    for mutr in ['5_prime_utr', '3_prime_utr', 'intronic_edge']:
                        val = stats[_file][stat][mutr][2]
                        try:
                            per = val * 100.0 / stats[_file][mutr][2]
                        except:
                            per = 0
                        lines[xx] += '{}({})\t'.format(val, round(per, rnd))
                        xx += 1

            if stat != 'SQNR':
                for mutr in ['5_prime_utr', '3_prime_utr', 'intronic_edge']:
                    val2 = stats[_file][stat][mutr][3]
                    try:
                        perc = len(val2) * 100.0 / len(stats[_file][mutr][3])
                    except:
                        perc = 0
                    try:
                        pert = sum(val2) * 100.0 / sum(stats[_file][mutr][3])
                    except:
                        pert = 0
                    sval2 = stats[_file][stat][mutr + '_u']
                    lines[xx] += '{}({})|{}({})\t'.format(sval2, round(perc, rnd),
                                                            sum(val2), round(pert, rnd))
                    xx += 1
                    mutr = mutr.replace('utr', 'ch').replace('edge', 'ch')
                    for k in range(2):
                        val = stats[_file][stat][mutr][k]
                        val = np.array(val)
                        vl = stats[_file][stat][mutr + '_u'][k]
                        try:
                            perc = vl * 100.0 / sval2
                        except:
                            perc = 0
                        try:
                            pert = sum(val) * 100.0 / sum(val2)
                        except:
                            pert = 0
                        lines[xx] += '{}({})|{}({})|{}|{}\t'.format(vl, round(perc, rnd),
                                                                    val.sum(), round(pert, rnd),
                                                                    round(val.mean(), rnd),
                                                                    round(val.std(), rnd))
                        xx += 1
                    val = stats[_file][stat][mutr + '_u'][2]
                    try:
                        perc = val * 100.0 / sval2
                    except:
                        perc = 0
                    lines[xx] += '{}({})\t'.format(val, round(perc, rnd))
                    xx += 1
            for MQR in ['MQ', 'MR']:
                if 'MQ' in stat and MQR == 'MQ':
                    gos = ['qgaps', 'qolaps']
                elif 'MR' in stat and MQR == 'MR':
                    gos = ['rgaps', 'rolaps']
                else:
                    continue
                for go in gos:
                    val2 = stats[_file][stat][go + 'c']
                    lines[xx] += '{}|{}\t'.format(len(val2), sum(val2))
                    xx += 1
                    for sstat in ['intronic', '3intron', 'intron5', '35prime']:
                        val = stats[_file][stat][go][sstat]
                        val = np.array(val)
                        if len(val2) != 0:
                            try:
                                perc = len(val) * 100.0 / len(val2)
                            except:
                                perc = 0
                            try:
                                pert = val.sum() * 100.0 / sum(val2)
                            except:
                                pert = 0
                        else:
                            perc = 0
                            pert = 0
                        lines[xx] += '{}({})|{}({})|{}|{}\t'.format(len(val), round(perc, rnd),
                                                                    val.sum(), round(pert, rnd),
                                                                    round(val.mean(), rnd),
                                                                    round(val.std(), rnd))
                        xx += 1
            # if 'MR' in stat:
            #     for go in ['rgaps', 'rolaps']:
            #         val2 = stats[_file][stat][go + 'c']
            #         lines[xx] += '{}\t'.format(val2)
            #         xx += 1
            #         for sstat in ['intronic', '3intron', 'intron5', '35prime']:
            #             val = stats[_file][stat][go][sstat]
            #             val = np.array(val)
            #             per = len(val) * 100.0 / val2
            #             lines[xx] += '{}({})|{}|{}|{}\t'.format(len(val), round(per, rnd), val.sum(),
            #                                                 round(val.mean(), rnd), round(val.std(), rnd))
            #             xx += 1
        lines[0] += '\t#\t'
        lines[1] += '\t#\t'
        lines[2] += '\t#\t'
        lines[3] += '\t#\t'

    for line in lines:
        write.write(line + '\n')

        # write.write('\n')
        # write.write('\n')
# val = stats['QUERY_TEST_SAMPLE.gtf.db.itracking']['count_u']
# print val
# val = np.array(val)
# print '{}|{}|{}|{}\t'.format(len(val), val.sum(),
#                                    round(val.mean(), 1), round(val.std(), 1))
# sm = 0
# for x in ['-', '+', '.']:
#     for b, e, d in stats['QUERY_TEST_SAMPLE.gtf.db.itracking']['ulength'][x]:
#         print b, e, e-b, x
#         sm += e-b
# print sm
