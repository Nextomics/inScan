#!/usr/bin/env python3
#author: archieyoung <yangqi2@grandomics.com>
import pysam
import sys


import cigar_parse


class fragment(object):
    def __init__(self,qname,ref,ref_start,strand,cigar,mapq):
        self.query_name = qname
        self.ref = ref
        self.ref_start = int(ref_start)
        self.strand = strand
        if self.strand == "+":
            self.cigar = cigar_parse.Cigar(cigar)
        else:
            self.cigar = cigar_parse.Cigar(cigar_parse.Cigar(cigar).reversed_cigar)
        self.mapq = mapq
        self.query_aligned_len = self.cigar.query_mapped_len
        self.ref_end = self.ref_start + self.cigar.ref_len - 1
        self.query_start = self.cigar.query_start
        self.query_end = self.cigar.query_end

class fragments(object):
    def __init__(self,sam_record):
        #sam is a pysam AlignmentFile
        self.sam_record = sam_record
        self.read_fragments = self._get_fragments()

    def _get_fragments(self):
        _fragments = []
        #'main' fragment
        qname = self.sam_record.query_name
        ref = self.sam_record.reference_name
        ref_start = self.sam_record.reference_start + 1
        mapq = self.sam_record.mapping_quality
        cigar = self.sam_record.cigarstring
        if self.sam_record.is_reverse:
            strand = "-"
        else:
            strand = "+"
        _fragment = fragment(qname,ref,ref_start,strand,cigar,mapq)
        _fragments.append(_fragment)
        #'other' fragment in SA:Z
        sa_tag = None
        try:
            sa_tag = self.sam_record.get_tag("SA")
        except KeyError:
            #not a split reads or sam format error
            pass
        if not sa_tag == None:
            sa_fields = sa_tag.split(";")[:-1]
            for i in sa_fields:
                ref, ref_start, strand, cigar, mapq = i.split(",")[:5]
                _fragment = fragment(qname,ref,ref_start,strand,cigar,mapq)
                _fragments.append(_fragment)
        return _fragments

def main():
    sam = pysam.AlignmentFile(sys.argv[1],"rb")
    for i in sam:
        read_fragments = fragments(i).read_fragments
        for j in read_fragments:
            print("\t".join([str(k) for k in [j.query_name,j.ref,j.ref_start,
                j.ref_end,j.query_start,j.query_end,j.strand,j.mapq,
                j.query_aligned_len]]))

if __name__ == "__main__":
    main()

