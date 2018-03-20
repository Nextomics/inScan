#!/usr/bin/env python3
#author archieyoung<yangqi2@grandomics.com>

"""
Find insertion from evidence of read split mapping.
Duplication is a special case of insertion, so I treat it as insertion here
"""

from itertools import combinations


class insertion(object):
    def __init__(self,query_name,ref,ref_start,ref_end,length):
        self.query_name = query_name
        self.ref = ref
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.length = length

class inter_ins(object):
    def __init__(self,read_fragments,min_len):
        self.read_fragments = sorted(read_fragments,key=lambda x:x.query_start)
        self.insertions = self._insertions()

    def _get_insert(self):
        if len(self.read_fragments) <= 1:
            return 0
        _insertions = []
        #A,B
        read_fragments_combine = combinations(self.read_fragments,2)
        for fr in read_fragments_combine:
            if fr[0].strand != fr[1].strand or fr[0].ref != fr[1].ref:
                return 0
            ins_length = ((fr[1].query_start-fr[0].query_end)-
                    (fr[1].ref_start-fr[0].ref_end))
            if ins_length >= min_len:
                if fr[1].ref_start > fr[0].ref_end:
                    ins_ref_start = fr[0].ref_end
                    ins_ref_end = fr[1].ref_start
                else:
                    ins_ref_start = fr[1].ref_start
                    ins_ref_end = fr[0].ref_end
            _ins = insertion(fr[0].query_name,fr[0].ref,ins_ref_start,
                    ins_ref_end,ins_length)
            _insertions.append(_ins)
        return _insertions
