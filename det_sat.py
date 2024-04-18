#!/usr/bin/env python

import sys
import os
import random
import argparse

#seqtk and mashmap should be in your path

def main():
    parser = argparse.ArgumentParser(description='Rough detection of _known_ satellite arrays, using previous annotation. seqtk and mashmap should be in your path.')
    parser.add_argument('--ref', help='Reference file', required=True)
    parser.add_argument('--bed', help='BED file with repeat annotation. Satelite subfamilies are considered to be the same; everything after first "_" or "(" in sat name will be ignored.',required=True)
    parser.add_argument('--asm', help='Assembly file', required=True)
    parser.add_argument('-o','--out_dir', help='Output directory', required=True)
    parser.add_argument('--chunk_len', type=int, default=10000, help='Chunk length')
    parser.add_argument('--fraction', type=float, default=0.05, help='Fraction of used repetitive chunks. For each array at least one chunk will be used.')
    parser.add_argument('--min_clust', type=int, default=100000, help='Minimal reported array size')
    parser.add_argument('--clust_dist', type=int, default=10000, help='Maximal allowed jump between mashmap alignments in a array')
    parser.add_argument('--jump_frac', type=float, default=0.2, help='Maximal allowed fraction of jumped length in a array')

    parser.add_argument('-t', type=int, default=20, help='threads')
    args = parser.parse_args()

    ref_file = args.ref
    bed_file = args.bed
    asm_file = args.asm
    out_dir = args.out_dir
    fraction = args.fraction
    chunk_len = args.chunk_len
    minclust = args.min_clust
    clust_dist = args.clust_dist
    jump_frac = args.jump_frac
    threads = args.t
    synonims={"dhor":"asat", "mon":"asat", "hor":"asat","LSAU-BSAT":"bsat"}
    synonims={}
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    random.seed(42)
    
    rep_bed = os.path.join(out_dir, "repsat.bed")
    rep_fasta = os.path.join(out_dir, "repsat.fasta")
    rep_tmp_fasta = os.path.join(out_dir, "repsat.fasta.tmp")

    res_mashmap = os.path.join(out_dir, 'mashmap.res')
    res_out = os.path.join(out_dir, 'detsat.out.bed')
    res_out_file = open(res_out, "w")
    rep_bed_file = open(rep_bed, "w")
    global_count = 0
    sat_names = {}
    for line in open(bed_file):
        arr = line.split()
        if len(arr) < 4:
            continue
        if not (arr[1].isdigit() and arr[2].isdigit()):
            continue
        #TODO: make it an option?
        sat_name = arr[3].split('_')[0].split('(')[0]
        sat_line = line.strip()
        sat_len = int(arr[2]) - int(arr[1])
   
        if sat_len < chunk_len:
            continue
#        print(f"Using array for {sat_name} with length {sat_len} bp")
        chunk_count = (int(arr[2]) - int(arr[1]))//chunk_len
        middle_pos = chunk_count//2
        array_start = int(arr[1])
        for i in range(0, chunk_count):
            if i == middle_pos or random.random() < fraction:
                chunk_start = i*chunk_len + array_start
                chunk_end = chunk_start + chunk_len
                if chunk_end > int(arr[2]):
                    chunk_end = int(arr[2])
                if chunk_end - chunk_start < chunk_len:
                    continue
                global_count += 1
                arr[1] = chunk_start
                arr[2] = chunk_end
                arr[3] += f"_chunk{global_count}"
                out_line = "\t".join(map(str, arr[:4])) + '\n'
                rep_bed_file.write(out_line)
                #seqtk specifics
                sat_names[f"{arr[0]}:{arr[1]+1}-{arr[2]}"] = sat_name
    rep_bed_file.close()
    os.system(f"rm {rep_fasta}")
    os.system(f"rm {rep_tmp_fasta}")
#    os.system (f"rm {res_mashmap}")
    #seen some problems with lowercase kmers, uppercasing everything
#    os.system(f"seqtk subseq {ref_file} {rep_bed} | tr '[:lower:]' '[:upper:]' > {rep_tmp_fasta}")
    os.system(f"seqtk subseq {ref_file} {rep_bed}  > {rep_tmp_fasta}")

    rep_fasta_file = open(rep_fasta, "w")
    for line in open(rep_tmp_fasta):
        if line[0] == '>':
            header = line.strip()[1:].split()[0]
            #print (header)
            line = ">" + sat_names[header]+"_" +line[1:].strip() + "\n"
        rep_fasta_file.write(line)
    rep_fasta_file.close()
    if os.path.exists(res_mashmap):
        print ("Reusing previous mashmap results")
    else:
        os.system(f"mashmap -r {rep_fasta} -q {asm_file} -o {res_mashmap} -s {chunk_len} -M -f map --dense --pi 85  -t {threads}")# 2> /dev/null")
    clusts = {}
    for line in open(res_mashmap):
        arr = line.split()
        if len(arr) < 12:
            continue
        #here we add processing of multiple sat representatives
        satname = arr[5].split('_')[0]


        if satname in synonims.keys():
            satname = synonims[satname]        
        contigname = arr[0]
        if contigname not in clusts:
            clusts[contigname] = {}
        if satname not in clusts[contigname]:
            clusts[contigname][satname] = []
        clusts[contigname][satname].append([int(arr[2]), int(arr[3])])    
    for contigs in sorted(clusts.keys()):
        lines = []
        for satname in sorted(clusts[contigs].keys()):
            clusts[contigs][satname].sort(key = lambda x: x[0])
            start_ind = 0
            while start_ind < len(clusts[contigs][satname]):
                end_ind = start_ind + 1
                jumped_len = 0
                start_pos = clusts[contigs][satname][start_ind][0]
                end_pos = clusts[contigs][satname][start_ind][1]
                if end_ind < len(clusts[contigs][satname]):
                    cur_jump = max(clusts[contigs][satname][end_ind][0] - end_pos, 0)                                        
                else:
                    cur_jump = 0
                #zero jumps always allowed, nonzeroes should be smaller than clust_dist and either total current length should be smaller than minclust or gap fraction less than jump_frac
                while end_ind < len(clusts[contigs][satname]) and ((cur_jump == 0) or (cur_jump <= clust_dist and ((jumped_len+ cur_jump)/(end_pos - start_pos) <= jump_frac or (end_pos - start_pos) < minclust))):
                    jumped_len += cur_jump
                    if end_pos < clusts[contigs][satname][end_ind][1]:
                        end_pos = clusts[contigs][satname][end_ind][1]
                    end_ind += 1
                    if end_ind < len(clusts[contigs][satname]):
                        cur_jump = max(clusts[contigs][satname][end_ind][0] - end_pos, 0)
                    else:
                        cur_jump = 0
                    
                if end_pos - start_pos >= minclust:
                    if jumped_len/(end_pos - start_pos) < jump_frac:
                        lines.append(f"{contigs}\t{start_pos}\t{end_pos}\t{satname}\t{end_pos - start_pos}\n")
                        print(f"Detected {satname} in {contigs} starting at {start_pos} of len {end_pos - start_pos} and {jumped_len/(end_pos - start_pos)}  fraction of jumped length")
                    else:
                        print(f"Detected {satname} in {contigs} starting at {start_pos} of len {end_pos - start_pos} and {jumped_len/(end_pos - start_pos)} fraction of jumped length, but skipped due to high jump fraction")
                start_ind = end_ind
        lines.sort(key = lambda x: int(x.split()[1]))
        for line in lines:
            res_out_file.write(line)
if __name__ == "__main__":
    main()
