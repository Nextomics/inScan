#!/usr/bin/env python3
#author:archieyoung<yangqi2@grandomics.com>
import sys
import json
import pysam


import fragment
import inter_ins
import intra_ins


def inScan(sam,min_ins):
    ins = {}
    #checked reads
    checked = []
    for record in sam:
        query_name = record.query_name
        if query_name not in checked:
            #read fragments
            read_fragments = fragment.fragments(record).read_fragments
            #inter fragments insertions
            inter_insertions = inter_ins.inter_ins(read_fragments,
                    min_ins).insertions
            if inter_insertions:
                ins.setdefault(query_name,
                        []).extend(inter_insertions)
        checked.append(query_name)
        #intra fragment insertions
        intra_insertions = intra_ins.intra_ins(record,min_ins).insertions
        if intra_insertions:
            ins.setdefault(query_name,[]).extend(intra_insertions)
    return ins

def overlap(a,b):
    if a[0] != b[0]:
        return 0
    if a[1] > b[2] or a[2] < b[1]:
        return 0
    else:
        return 1

def regions_inScan(sam_io,bed,min_ins):
    region_ins_all = {}
    for region in bed:
        #fetch sam
        sam = sam_io.fetch(region[0],region[1],region[2])
        ins = inScan(sam,min_ins)
        ins_in_region = {}
        for query_name in ins:
            keep = []
            for i in ins[query_name]:
                if overlap([i.ref,i.ref_start,i.ref_end],
                        [region[0],region[1],region[2]]):
                    #print(i.query_name,i.ref,i.ref_start,i.ref_end,i.length)
                    keep.append([i.ref,i.ref_start,i.ref_end,i.length])
            if keep == []:
                continue
            ins_in_region.setdefault(query_name,keep)
            region_ins_all.setdefault(region[0]+":"+str(region[1])+"-"+str(region[2])
                    ,ins_in_region)
    return region_ins_all

def main():
    sam_io = pysam.AlignmentFile(sys.argv[1],"rb")
    bed = []
    with open(sys.argv[2],"r") as bed_io:
        for line in bed_io:
            fields = line.strip().split("\t")
            bed.append((fields[0],int(fields[1]),int(fields[2])))
    region_ins_all = regions_inScan(sam_io,bed,20)
    with open(sys.argv[3],"w") as out_json:
        json.dump(region_ins_all, out_json, sort_keys = False, indent=4)

if __name__ == "__main__":
    main()
