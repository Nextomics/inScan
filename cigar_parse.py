#!/usr/bin/env python3
#author: yangqi
#contact: yangqi2@grandomics.com
"""
cigar string is used to represent sequence alignment. It is compact.
cigar is used in sam format file and other files.
For more information about cigar or sam files see: https://samtools.github.io/hts-specs

Op BAM                 Description                             Consumes_query Consumes_reference
M   0  alignment match (can be a sequence match or mismatch)        yes            yes
I   1  insertion to the reference                                   yes             no
D   2  deletion from the reference                                   no            yes
N   3  skipped region from the reference                             no            yes
S   4  soft clipping (clipped sequences present in SEQ)             yes             no
H   5  hard clipping (clipped sequences NOT present in SEQ)          no             no
P   6  padding (silent deletion from padded reference)               no             no
=   7  sequence match                                               yes            yes
X   8  sequence mismatch                                            yes            yes

aligned_len
ref    ATTTGC-CGA
query  ATT-GCCCGA
"""
#qstart and qend need to be modified, 'fisrt M' as start or end pos

class Cigar(object):

    def __init__(self, cigar_string):
        self.CIGAR = cigar_string
        #consums mean step right ->
        self.query_consumes = ["M","I","S","=","X"]
        self.ref_consumes = ["M","D","N","=","X"]
        num = ""
        op = ""
        self.cigar_ops = []
        for c in self.CIGAR:
            if c.isdigit():
                num += c
            else:
                op = c
                self.cigar_ops.append((int(num), op))
                num = ""
                op = ""
        self.query_len = self._query_len()
        self.query_mapped_len = self._query_mapped_len()
        self.ref_len = self._ref_len()
        self.aligned_len = self._aligned_len()
        self.query_start = self._query_start()
        self.query_end = self._query_end()
        self.reversed_cigar = self._cigar_reverse()

    def _query_len(self):
        #query_len = num of hard-clip bases + num of query_consumes bases
        query_len_operations = self.query_consumes + ["H"]
        return sum([co[0] for co in self.cigar_ops \
                if co[1] in query_len_operations])

    def _query_mapped_len(self):
        query_mapped_len_operations = ["M","I","=","X"]
        return sum([co[0] for co in self.cigar_ops \
                if co[1] in query_mapped_len_operations])

    def _ref_len(self):
        return sum([co[0] for co in self.cigar_ops \
                if co[1] in self.ref_consumes])

    def _aligned_len(self):
        #query_mapped_len = except "H"s and "S"s
        not_aligned_len_operations = ["H", "S"]
        return sum([co[0] for co in self.cigar_ops \
                if co[1] not in not_aligned_len_operations])

    def _query_start(self):
        #one based start
        #"H" and "S" appear togather is not considered, eg "10H5S89M"
        #suit for "10H89M" or "10S89M"
        if self.cigar_ops[0][1] == "H" or self.cigar_ops[0][1] == "S":
            return 1+self.cigar_ops[0][0]
        else:
            return 1

    def _query_end(self):
        #one based end
        #"H" and "S" appear togather is not considered, eg "89M5S10H"
        #suit for "89M10H" or "89M10S"
        if self.cigar_ops[-1][1] == "H" or self.cigar_ops[-1][1] == "S":
            return self._query_len()-self.cigar_ops[-1][0]
        else:
            return self._query_len()

    def _cigar_reverse(self):
        return "".join(["".join([str(co[0]),co[1]]) for co in self.cigar_ops][::-1])


