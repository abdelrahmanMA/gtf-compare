#!/usr/bin/env python
# -*- coding: utf-8 -*-

class GtfInterval:
    '''
    The GtfInterval object represents a GTF Interval

    Parameters
    ----------
    interval : Interval
        This is the Interval object for this GTF interval

    exon :
        This is a reference to Exon Object

    number:
        This is the exon's order inside the parent transcript

    Attributes
    ----------
    interval : Interval
        This is a reference to this GTF Interval

    begin : int
        This the interval's start position

    end : int
        This the interval's end position

    strand : str
        This is the interval's strand

    chrom : str
        This is the interval's parent chromosome name

    geneIds : {str}
        This set contains all possible gene ids that contain this interval

    transcriptIds : {str: int}
        This dictionary contains all possible transcript ids that contain this interval with the order of the interval in the transcript

    exonIds : {str}
        This set contains all possible exon ids of this interval

    note : number
        This is the number of the matching note in predefined set of notes, to be reported later when generating files
    '''

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
        '''Adds all possible exon ids and transcript ids and gene ids

        Parameters
        ----------
        exon_id : str
            This is an alternative exon id for the interval

        transcript_id : str
            This is an alternative parent transcript id for the interval

        number : int
            This order of the interval in the alternative parent transcript

        gene_id : str
            This is an alternative parent gene id for the interval

        Raises
        ------
        RuntimeError
            Duplicate Transcript

        Returns
        -------
        None
        '''
        if exon_id not in self.exonIds:
            self.exonIds.add(exon_id)
        else:
            raise RuntimeError('Duplicate feature, Exon')

        if transcript_id not in self.transcriptIds:
            self.transcriptIds[transcript_id] = number
        else:
            raise RuntimeError('Duplicate Transcript')

        if gene_id not in self.geneIds:
            self.geneIds.add(gene_id)


class GtfExon:
    '''
    The GtfInterval object represents a GTF Exon

    Parameters
    ----------
    _id : str
        This is the exon's id

    gtf_interval : GtfInterval
        This is a reference to Gtf Interval instance

    transcript_id: str
        This is the exon's parent transcript id

    gene_id: str
        This is the exon's parent gene id

    Attributes
    ----------
    id : str
        This is the exon's id

    begin : int
        This the exon's start position

    end : int
        This the exon's end position

    strand : str
        This is the exon's strand

    chrom : str
        This is the exon's parent chromosome name

   gtf_interval : GtfInterval
        This is a reference to Gtf interval

    transcript_id: str
        This is the exon's parent transcript id

    gene_id: str
        This is the exon's parent gene id
    '''

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
    '''
    The GtfInterval object represents a GTF Transcript

    Parameters
    ----------
    _id : str
        This is the transcript's id

    gtf_exon : GtfExon
        This is a reference to Gtf Exon instance

    gene_id: str
        This is the transcripts's parent gene id

    Attributes
    ----------
    id : str
        This is the transcript's id

    begin : int
        This the transcript's start position

    end : int
        This the transcript's end position

   gtf_exons : {GtfExon}
        This set holds references to all Gtf Exon instances found in the transcript

    exonIds: {str}
        This set holds all exon ids found in the transcript

    gene_id: str
        This is the transcript's parent gene id
    '''

    def __init__(self, _id, gtf_exon, gene_id):
        interval = gtf_exon.gtf_interval.interval
        self.id = _id
        self.gtf_exons = {gtf_exon}
        self.exonIds = {gtf_exon.id}
        self.gene_id = gene_id
        self.begin = interval.begin
        self.end = interval.end

    def add(self, gtf_exon):
        '''Adds Exons to GTF Transcript

        Parameters
        ----------
        gtf_exon : GtfExon
            This is a reference to Gtf Exon instance to be added

        Returns
        -------
        None
        '''
        interval = gtf_exon.gtf_interval.interval
        self.gtf_exons.add(gtf_exon)
        self.exonIds.add(gtf_exon.id)
        self.begin = min(self.begin, interval.begin)
        self.end = max(self.end, interval.end)


class GtfGene:
    '''
    The GtfGene object represents a GTF Gene

    Parameters
    ----------
    _id : str
        This is the gene's id

    gtf_transcript: str
        This a reference to GTF Transcript instance

    Attributes
    ----------
    id : str
        This is the gene's id

    begin : int
        This the gene's start position

    end : int
        This the gene's end position

    gtf_transcripts : {GtfTranscript}
        This set holds references to all Gtf Transcript instances found in the gene

    transcriptsIds: {str}
        This set holds all transcript ids found in the gene

    exonIds: {str}
        This set holds all exon ids found in the gene

    '''

    def __init__(self, _id, gtf_transcript):
        self.id = _id
        self.gtf_transcripts = {gtf_transcript}
        self.transcriptIds = {gtf_transcript.id}
        self.exonIds = gtf_transcript.exonIds
        self.begin = gtf_transcript.begin
        self.end = gtf_transcript.end

    def add(self, gtf_transcript):
        '''Adds Transcripts to GTF Gene

        Parameters
        ----------
        gtf_transcript : GtfTranscript
            This is a reference to Gtf Transcript instance to be added

        Returns
        -------
        None
        '''

        self.gtf_transcripts.add(gtf_transcript)
        self.transcriptIds.add(gtf_transcript.id)
        self.exonIds |= gtf_transcript.exonIds
        self.begin = min(self.begin, gtf_transcript.begin)
        self.end = max(self.end, gtf_transcript.end)
