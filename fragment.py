#!/usr/bin/env
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
    def __init__(self,sam_record,sam):
        #sam is a pysam AlignmentFile
        self.sam = sam
        self.sam_record = sam_record
        self.read_fragments = self._get_fragments()

    def _get_fragments(self):
        _fragments = []
        sam_str = self.sam_record.tostring(self.sam)
        #'main' fragment
        fields = sam_str.strip().split("\t")
        (qname, flag, ref, ref_start, mapq,
                cigar, mrnm, mpos, tlen, seq, qual) = fields[:11]
        if self.sam_record.is_reverse:
            strand = "-"
        else:
            strand = "+"
        _fragment = fragment(qname,ref,ref_start,strand,cigar,mapq)
        _fragments.append(_fragment)
        #optional tags
        optional_tags = fields[11:]
        #'other' fragment in SA:Z
        sa_tag = None
        for tag in optional_tags:
            if tag[:3] == "SA:Z":
                sa_tag = tag[3:]
        if not sa_tag == None:
            sa_fields = sa_tag.split(";")
            for i in sa_fields:
                ref, ref_start, strand, cigar, mapq = i.split(",")
                _fragment = fragment(qname,ref,ref_start,strand,cigar,mapq)
                _fragments.append(_fragment)
        return _fragments

def main():
    sam = pysam.AlignmentFile(sys.argv[1],"rb")
    for i in sam:
        read_fragments = fragments(i,sam).read_fragments
        for j in read_fragments:
            print("\t".join([str(k) for k in [j.query_name,j.ref,j.ref_start,
                j.ref_end,j.query_start,j.query_end,j.strand,j.mapq,
                j.query_aligned_len]]))

if __name__ == "__main__":
    main()

