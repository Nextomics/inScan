#!/usr/bin/env python3
#author: archieyoung<yangqi2@grandomics.com>
"""
find insertions from cigar
"""

import cigar_parse
from inter_ins import insertion

class intra_ins(object):
    def __init__(self,sam_record):
