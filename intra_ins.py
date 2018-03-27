#!/usr/bin/env python3
#author: archieyoung<yangqi2@grandomics.com>
"""
find insertions from cigar
"""

import cigar_parse
from inter_ins import insertion


class intra_ins(object):
    def __init__(self,sam_record,min_ins_len):
        self.query_name = sam_record.query_name
        #minimum sv length to take into account
        self.min_ins_len = min_ins_len
        #change to 1 based coordinate
        self.reference_start = sam_record.reference_start
        self.reference_name = sam_record.reference_name
        #cigar_parse instance
        self.cigar = cigar_parse.Cigar(sam_record.cigarstring)
        self.insertions = self._get_insert()

    def _get_insert(self):
        _insertions = []
        ref_pos = self.reference_start
        for num,op in self.cigar.cigar_ops:
            if op in self.cigar.ref_consumes:
                ref_pos += num
            if op == "I" and num >= self.min_ins_len:
                _ins = insertion(self.query_name,self.reference_name,
                        ref_pos,ref_pos,num)
                _insertions.append(_ins)
        return _insertions

