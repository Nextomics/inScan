#!/usr/bin/env
#author: archieyoung <yangqi2@grandomics.com>

class fragment(object):
    def __init__(self,ref,ref_start,strand,cigar,mapq,query_aligned_len):
        _fields = sam_record.split("\t")
        self.ref = ref
        self.ref_start = int(ref_start)
        self.strand = strand
        self.mapq = mapq
        self.query_aligned_len = query_aligned_len

        self.ref_end =

