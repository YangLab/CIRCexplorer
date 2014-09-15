#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
CIRCexplorer.py 1.0.0 -- circular RNA analysis toolkits.

Usage: CIRCexplorer.py [options]

Options:
    -f FUSION, --fusion=FUSION      Fusion BAM file
    -g GENOME, --genome=GENOME      Genome FASTA file
    -r REF, --ref=REF               Gene annotation
    -o OUT, --output=OUT            Output file [default: circ.txt]
"""

__author__ = 'Xiao-Ou Zhang (zhangxiaoou@picb.ac.cn)'
__version__ = '1.0.0'

from docopt import docopt
import sys
import pysam
from collections import defaultdict
from interval import Interval
import tempfile
import os

def convert_fusion(fusion_bam, output):
    """
    Extract fusion junction reads from the BAM file
    """
    print('Start to convert fustion reads...')
    fusion = defaultdict(int)
    for i, read in enumerate(parse_bam(fusion_bam)):
        chrom, strand, start, end = read
        segments = [start, end]
        if (i + 1) % 2 == 1:  # first fragment of the fusion junction read
            interval = [start, end]
        else:  # second fragment of the fusion junction read
            sta1 = interval[0]
            end1 = interval[1]
            sta2 = segments[0]
            end2 = segments[1]
            if end1 < sta2 or end2 < sta1:  # no overlap between fragments
                sta = sta1 if sta1 < sta2 else sta2
                end = end1 if end1 > end2 else end2
                fusion['%s\t%d\t%d' % (chrom, sta, end)] += 1
    total = 0
    with open(output, 'w') as outf:
        for i, pos in enumerate(fusion):
            outf.write('%s\tFUSIONJUNC_%d/%d\t0\t+\n' % (pos, i, fusion[pos]))
            total += fusion[pos]
    print('Converted %d fusion reads!' % total)

def annotate_fusion(ref_f, input, output):
    """
    Align fusion juncrions to gene annotations
    """
    print('Start to annotate fusion junctions...')
    gene, isoform = parse_ref1(ref_f) # gene annotations
    fusion, fusion_index = parse_bed(input) # fusion junctions
    total = 0
    with open(output, 'w') as outf:
        for chrom in gene:
            # overlap gene annotations with fusion juncrions
            result = Interval.overlapwith(gene[chrom].interval, fusion[chrom])
            for itl in result:
                # extract isoform annotations
                iso = list(filter(lambda x: x.startswith('iso'), itl[2:]))
                for fus in itl[(2 + len(iso)):]: # for each overlapped fusion junction
                    mapper, name, reads = fus.split()[1:]
                    fus_start, fus_end = fusion_index[fus]
                    edge_annotation = '' # first or last exon flag
                    for iso_id in iso:
                        g, i, c, s = iso_id.split()[1:]
                        start = isoform[iso_id][0][0]
                        end = isoform[iso_id][-1][-1]
                        if fus_start < start or fus_end > end: # fusion junction excesses boundaries of isoform annotation
                            continue
                        fusion_info, index, edge = map_fusion_to_iso(fus_start,
                                                                     fus_end, s,
                                                                     isoform[iso_id])
                        if fusion_info:
                            fus_start_str = str(fus_start)
                            fus_end_str = str(fus_end)
                            bed_info = '\t'.join([chrom, fus_start_str,
                                                  fus_end_str,
                                                  '/'.join([mapper, reads]),
                                                  '0', s, fus_start_str,
                                                  fus_start_str, '0,0,0'])
                            bed = '\t'.join([bed_info, fusion_info, g, i, index])
                            if not edge: # not first or last exon
                                outf.write(bed + '\n')
                                total += 1
                                break
                            else: # first or last exon
                                edge_annotation = bed
                    else: # cannot align or boundary exon
                        if edge_annotation: # first or last exon
                            outf.write(bed + '\n')
                            total += 1
    print('Annotated %d fusion junctions!' % total)

def fix_fusion(ref_f, genome_fa, input, output):
    """
    Realign fusion juncrions
    """
    print('Start to fix fusion junctions...')
    fa = genome_fa
    ref = parse_ref2(ref_f)
    fusion, fixed_flag = fix_bed(input, ref, fa)
    total = 0
    with open(output, 'w') as outf:
        for fus in fusion:
            tophat_reads = fusion[fus]['tophat']
            fixed = str(fixed_flag[fus]['tophat'])
            if fixed == '1':
                total += 1
            name = 'circular_RNA/' + str(tophat_reads)
            gene, iso, chrom, strand, index = fus.split()
            starts, ends = ref['\t'.join([gene, iso, chrom, strand])]
            if ',' in index:  # back spliced exons
                s, e = [int(x) for x in index.split(',')]
                start = str(starts[s])
                end = str(ends[e])
                length = str(e - s + 1)
                sizes, offsets = generate_bed(int(start), starts[s:(e + 1)],
                                              ends[s:(e + 1)])
                if s == 0:
                    left_intron = 'None'
                else:
                    left_intron = '%s:%d-%d' % (chrom, ends[s - 1], starts[s])
                if e == len(ends) - 1:
                    right_intron = 'None'
                else:
                    right_intron = '%s:%d-%d' % (chrom, ends[e], starts[e + 1])
                intron = '|'.join([left_intron, right_intron])
                bed = '\t'.join([chrom, start, end, name, fixed, strand, start,
                                 start, '0,0,0', length, sizes, offsets,
                                 str(tophat_reads), 'No', gene, iso, intron])
            else:  # ciRNAs
                index, start, end = index.split('|')
                size = str(int(end) - int(start))
                index = int(index)
                intron = '%s:%d-%d' % (chrom, ends[index], starts[index + 1])
                bed = '\t'.join([chrom, start, end, name, fixed, strand, start,
                                 start, '0,0,0', '1', size, '0',
                                 str(tophat_reads), 'Yes', gene, iso, intron])
            outf.write(bed + '\n')
    print('Fixed %d fusion junctions!' % total)

def parse_bam(bam):
    fusion = {}
    for read in bam:
        if read.is_secondary: # not the primary alignment
            continue
        tags = dict(read.tags)
        if 'XF' not in tags: # not fusion junctions
            continue
        chr1, chr2 = tags['XF'].split()[1].split('-')
        if chr1 != chr2: # not on the same chromosome
            continue
        strand = '+' if not read.is_reverse else '-'
        if read.qname not in fusion: # first fragment
            fusion[read.qname] = [chr1, strand, read.pos, read.aend]
        else: # second fragment
            if chr1 == fusion[read.qname][0] and strand == fusion[read.qname][1]:
                yield [chr1, strand, read.pos, read.aend]
                yield fusion[read.qname]

def parse_ref1(ref_file):
    gene = defaultdict(list)
    iso = {}
    with open(ref_file, 'r') as f:
        for line in f:
            gene_id, iso_id, chrom, strand = line.split()[:4]
            total_id = '\t'.join(['iso', gene_id, iso_id, chrom, strand])
            starts = list(map(int, line.split()[9].split(',')[:-1]))
            ends = list(map(int, line.split()[10].split(',')[:-1]))
            start = starts[0]
            end = ends[-1]
            gene[chrom].append([start, end, total_id])
            iso[total_id] = [starts, ends]
    for chrom in gene:
        gene[chrom] = Interval(gene[chrom])
    return (gene, iso)

def parse_ref2(ref_file):
    gene = {}
    with open(ref_file, 'r') as f:
        for line in f:
            gene_id, iso_id, chrom, strand = line.split()[:4]
            starts = [int(x) for x in line.split()[9].split(',')[:-1]]
            ends = [int(x) for x in line.split()[10].split(',')[:-1]]
            gene['\t'.join([gene_id, iso_id, chrom, strand])] = [starts, ends]
    return gene

def parse_bed(fus):
    fusion = defaultdict(list)
    fusion_index = {}
    with open(fus, 'r') as f:
        for line in f:
            chrom, start, end, name = line.split()[:4]
            start = int(start)
            end = int(end)
            reads = name.split('/')[1]
            fusion_id = 'fusion\ttophat\t%s\t%s' % (name, reads)
            fusion[chrom].append([start, end, fusion_id])
            fusion_index[fusion_id] = [start, end]
    return (fusion, fusion_index)

def map_fusion_to_iso(start, end, strand, iso_info):
    starts = iso_info[0]
    ends = iso_info[1]
    # check sequnence within +/-10bp
    start_points = list(range(start - 10, start + 11))
    end_points = list(range(end - 10, end + 11))
    start_index, end_index = None, None
    start_intron_flag, end_intron_flag = False, False
    # check starts
    for i, s in enumerate(starts):
        if s in start_points:
            start_index = i
            break
    else:
        for j, e in enumerate(ends):
            if e in start_points:
                start_index = j
                start_intron_flag = True
                break
    # check ends
    for j, e in enumerate(ends):
        if e in end_points:
            end_index = j
            break
    else:
        for i, s in enumerate(starts):
            if s in end_points:
                end_index = i - 1
                end_intron_flag = True
                break
    # ciRNAs
    if start_intron_flag and strand == '+' and end < starts[start_index + 1]:
        return ('\t'.join(['1', str(end - start), '0', 'Yes']),
                str(start_index), False)
    elif end_intron_flag and strand == '-' and start > ends[end_index]:
        return ('\t'.join(['1', str(end - start), '0', 'Yes']), str(end_index),
                False)
    # back spliced exons
    elif (start_index is not None and end_index is not None and
          not start_intron_flag and not end_intron_flag):
        if start_index != 0 and end_index != len(ends) - 1:
            edge_flag = False
        else:
            edge_flag = True
        return convert_to_bed(start, end,
                            starts[start_index:(end_index + 1)],
                            ends[start_index:(end_index + 1)],
                            str(start_index), str(end_index), edge_flag)
    else:
        return(None, None, False)

def convert_to_bed(start, end, starts, ends, start_index, end_index, edge):
    new_starts = [start] + starts[1:]
    new_ends = ends[:-1] + [end]
    block_starts, block_sizes = [], []
    for s, e in zip(new_starts, new_ends):
        block_starts.append(str(s - start))
        block_sizes.append(str(e - s))
    length = len(block_sizes)
    block_starts = ','.join(block_starts)
    block_sizes = ','.join(block_sizes)
    index = ','.join([start_index, end_index])
    return ('\t'.join([str(length), block_sizes, block_starts, 'No']), index,
            edge)

def fix_bed(fusion_file, ref, fa):
    fusion = defaultdict(dict)
    fixed_flag = defaultdict(dict) # flag to indicate realignment
    with open(fusion_file, 'r') as f:
        for line in f:
            chrom = line.split()[0]
            strand = line.split()[5]
            start, end = [int(x) for x in line.split()[1:3]]
            mapper, reads = line.split()[3].split('/')
            reads = int(reads)
            flag, gene, iso, index = line.split()[-4:]
            flag = True if flag == 'Yes' else False
            name = '\t'.join([gene, iso, chrom, strand, index])
            iso_starts, iso_ends = ref['\t'.join([gene, iso, chrom, strand])]
            if not flag: # back spliced exons
                s, e = [int(x) for x in index.split(',')]
                if start == iso_starts[s] and end == iso_ends[e]: # not realign
                    fusion[name][mapper] = fusion[name].get(mapper, 0) + reads
                    fixed_flag[name][mapper] = 0
                elif check_seq(chrom, [start, iso_starts[s], end, iso_ends[e]],
                               fa): # realign
                    fusion[name][mapper] = fusion[name].get(mapper, 0) + reads
                    fixed_flag[name][mapper] = 1
            else: # ciRNAs
                index = int(index)
                if strand == '+':
                    if start == iso_ends[index]: # not realign
                        name += '|'.join(['', str(start), str(end)])
                        fusion[name][mapper] = fusion[name].get(mapper, 0) + reads
                        fixed_flag[name][mapper] = 0
                    elif check_seq(chrom, [start, iso_ends[index], end], fa,
                                   intron_flag=True): # realign
                        fixed_start = iso_ends[index]
                        fixed_end = end + fixed_start - start
                        name += '|'.join(['', str(fixed_start), str(fixed_end)])
                        fusion[name][mapper] = fusion[name].get(mapper, 0) + reads
                        fixed_flag[name][mapper] = 1
                else:
                    if end == iso_starts[index + 1]: # not realign
                        name += '|'.join(['', str(start), str(end)])
                        fusion[name][mapper] = fusion[name].get(mapper, 0) + reads
                        fixed_flag[name][mapper] = 0
                    elif check_seq(chrom, [end, iso_starts[index + 1], start],
                                   fa, intron_flag=True): # realign
                        fixed_end = iso_starts[index + 1]
                        fixed_start = start + fixed_end - end
                        name += '|'.join(['', str(fixed_start), str(fixed_end)])
                        fusion[name][mapper] = fusion[name].get(mapper, 0) + reads
                        fixed_flag[name][mapper] = 1
    return (fusion, fixed_flag)

def check_seq(chrom, pos, fa, intron_flag=False):
    if not intron_flag: # back spliced exons
        if pos[0] - pos[1] != pos[2] - pos[3]:
            return False
        if pos[0] < pos[1]:
            seq1 = fa.fetch(chrom, pos[0], pos[1])
            seq2 = fa.fetch(chrom, pos[2], pos[3])
        else:
            seq1 = fa.fetch(chrom, pos[1], pos[0])
            seq2 = fa.fetch(chrom, pos[3], pos[2])
    else: # ciRNAs
        if abs(pos[0] - pos[1]) <= 5:  # permit mismatches within 5bp
            return True
        elif pos[0] < pos[1]:
            seq1 = fa.fetch(chrom, pos[0], pos[1])
            seq2 = fa.fetch(chrom, pos[2], pos[2] + pos[1] - pos[0])
        else:
            seq1 = fa.fetch(chrom, pos[1], pos[0])
            seq2 = fa.fetch(chrom, pos[2] - pos[0] + pos[1], pos[2])
    if seq1 == seq2:
        return True
    else:
        return False

def generate_bed(start, starts, ends):
    sizes, offsets = [], []
    for s, e in zip(starts, ends):
        sizes.append(str(e - s))
        offsets.append(str(s - start))
    sizes = ','.join(sizes)
    offsets = ','.join(offsets)
    return (sizes, offsets)

def create_temp():
    temp_dir = tempfile.mkdtemp()
    temp1 = temp_dir + '/tmp1'
    temp2 = temp_dir + '/tmp2'
    return (temp_dir, temp1, temp2)

def delete_temp(temp_dir, temp1, temp2):
    os.remove(temp1)
    os.remove(temp2)
    os.rmdir(temp_dir)

if __name__ == '__main__':
    options = docopt(__doc__)
    for arg in options:
        if not options[arg]:
            sys.exit(__doc__)
    try:
        fusion_bam = pysam.Samfile(options['--fusion'], 'rb')
    except:
        sys.exit('Please make sure %s is a BAM file!' % options['--fusion'])
    try:
        genome_fa = pysam.Fastafile(options['--genome'])
    except:
        sys.exit('Please make sure %s is a Fasta file and indexed!' % options['--genome'])
    ref_f = options['--ref']
    output = options['--output']
    temp_dir, temp1, temp2 = create_temp()
    convert_fusion(fusion_bam, temp1)
    annotate_fusion(ref_f, temp1, temp2)
    fix_fusion(ref_f, genome_fa, temp2, output)
    delete_temp(temp_dir, temp1, temp2)
