"""
author: archieyoung<yangqi2@grandomics.com>

Find insertion from evidence of read split mapping.
Duplication is a special case of insertion, so I treat it as insertion here
"""


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
        self.min_len = min_len
        self.insertions = self._get_insert()

    def _has_ins(self,fr1,fr2):
        if fr1.ref != fr2.ref or fr1.strand != fr2.strand:
            return 0
        #then fr1 and fr2 are in the same ref and hava the same strand
        if fr1.strand == "+":
            #fr_relativa_dist = fragment dist on reads - fragment dist on ref
            fr_relativa_dist = ((fr2.query_start-fr1.query_end)-
                    (fr2.ref_start-fr1.ref_end))
        else:
            fr_relativa_dist = ((fr2.query_start-fr1.query_end)-
                    (fr1.ref_start-fr2.ref_end))
        if fr_relativa_dist > 0:
            #has insertion
            return fr_relativa_dist
        else:
            return 0

    def _is_same_ref(self,fr1,fr2):
        #fragment1 and fragment2 are in the same ref
        if fr1.ref == fr2.ref:
            return 1
        else:
            return 0

    def _calculate_ins(self,fr1,fr2,fr_relativa_dist):
        if fr1.strand == "+":
            #define ins_len and asign 0 to it
            ins_len = 0
            #insertion
            #case1 and case2
            if fr2.ref_start >= fr1.ref_end:
                ins_len = fr2.query_start - fr1.query_end
                ins_ref_start = fr1.ref_end
                ins_ref_end = fr2.ref_start
            #case3
            if (fr2.ref_start < fr1.ref_end and
                    fr2.ref_end >= fr1.ref_end):
                ins_len = fr_relativa_dist
                ins_ref_start = fr1.ref_end
                ins_ref_end = fr1.ref_end
            #case4
            if (fr2.ref_start < fr1.ref_end and
                    fr2.ref_end < fr1.ref_end):
                ins_len = fr2.query_end - fr1.query_end
                ins_ref_start = fr1.ref_end
                ins_ref_end = fr1.ref_end
        else:
            #strand are "-"
            #define ins_len and asign 0 to it
            ins_len = 0
            #insertion
            #case1 or case2
            if fr1.ref_start >= fr2.ref_end:
                ins_len = fr2.query_start - fr1.query_end
                ins_ref_start = fr2.ref_end
                ins_ref_end = fr1.ref_start
            #case3
            if (fr1.ref_start < fr2.ref_end and
                    fr1.ref_end >= fr2.ref_end):
                ins_len = fr_relativa_dist
                ins_ref_start = fr2.ref_end
                ins_ref_end = fr2.ref_end
            #case4
            if (fr1.ref_start < fr2.ref_end and
                    fr1.ref_end < fr2.ref_end):
                ins_len = fr2.query_end - fr1.query_end
                ins_ref_start = fr2.ref_end
                ins_ref_end = fr2.ref_end

        _ins = insertion(fr1.query_name,fr1.ref,ins_ref_start,
                ins_ref_end,ins_len)

        return _ins


    def _get_insert(self):
        if len(self.read_fragments) <= 1:
            return 0

        _insertions = []

        #a new read fragment combination strategy, fix bug of false positive
        #insertion discovery
        while self.read_fragments:
            for fr in self.read_fragments[1:]:
                has_ins = self._has_ins(self.read_fragments[0],fr)
                if has_ins:
                    _ins = self._calculate_ins(self.read_fragments[0],
                            fr,has_ins)
                    if _ins.length >= self.min_len:
                        _insertions.append(_ins)
                        break
                else:
                    if self._is_same_ref(self.read_fragments[0],fr):
                        break
            self.read_fragments.pop(0)

        return _insertions
