from argparse import ArgumentParser

def argument_parser():
    '''Parses Agruments from Terminal'''
    parser = ArgumentParser()

    parser.add_argument('-r', dest='refgtf', help='''“reference” annotation GTF file. Each sample is matched against this file, and sample isoforms are tagged as overlapping, matching, or novel where appropriate. <input.gtf>''',
                        metavar='<reference_mrna>', nargs='?', required=True)

    parser.add_argument('-i', help='''provide a text file with a list of GTF files to process instead
       of expecting them as command line arguments (useful when a large number
       of GTF files should be processed)''', metavar='<input_list>', nargs='?')

    parser.add_argument('input_annotations', help='''input GTF files to be matched against the reference GTF<input.gtf>''', metavar='<input.gtf>', nargs='*')

    # parser.add_argument('-th', help='a basepair threshold for exon matching',
    #                     metavar='basepair threshold', type=int, nargs=1, default=[0])

    args = parser.parse_args()
    if len(args.input_annotations) == 0 and args.i is None:
        print('Not enough arguments, enter at least one assembled annotation file')

    input_annotations = []

    for x in args.input_annotations:
        input_annotations.append(x)
    if args.i is not None:
        with open(args.i, 'r') as input_list:
            for line in input_list:
                input_annotations.append(line.rstrip())

    input_files = input_annotations
    reference = str(args.refgtf)
    # threshold = args.th[0]
    return [reference, input_files]
