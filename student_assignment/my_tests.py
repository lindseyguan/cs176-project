""" 
    RNA Alignment Assignment

    My own tests
"""

import sys # DO NOT EDIT THIS
import ctypes
import ctypes.util
import functools
import numpy as np
import time
from tqdm import tqdm
from shared import *
from evaluation import *
from project import *

def testBWTFunctions():
    # s = 'ACTGGTTACCCTACTGATTAGGACTC$'
    s = 'ACACT$'
    # print(s)
    sa = get_suffix_array(s)
    for i in range(len(sa)):
        print(str(i) + ': ' + str(sa[i]) + ' -> ' + s[sa[i]:])
    L = get_bwt(s, sa)
    # print(L)
    F = get_F(L)
    # print(F)
    M = get_M(F)
    # print(M)
    occ = get_occ(L)
    # print(occ)
    matches = exact_suffix_matches('G', M, occ)
    print(matches)

def testAlignerInit():
    # Testing runtime of Aligner init
    # genome_sequence = ''
    # with open('./genome.fa') as f:
    #     genome_sequence = f.readline()
    #     genome_sequence = f.readline() + '$'
    genome_sequence = 'ACTGGTTACCCTACTGACCG'
    # read = 'ATTACTCTTGGGAATGAAATCCTATCTATATAAGCTGTGGTTTGAAATCC'
    read = 'ACTGCACCC'

    genes = set()

    gene_id = ''
    isoforms = []

    isoform_id = ''
    exons = []

    exon_id = ''
    start = 0
    end = 0
        
    for line in reversed(list(open("./genes_short.tab"))):
        elements = line.split('\t')
        if elements[0] == 'exon':
            ex = Exon(elements[1], int(elements[2]), int(elements[3]))
            exons.append(ex)
        elif elements[0] == 'isoform':
            iso = Isoform(elements[1], exons)
            isoforms.append(iso)
            exons = []
        elif elements[0] == 'gene':
            g = Gene(elements[1], isoforms)
            genes.add(g)
            isoforms = []

    # reads = []
    # for line in list(open("./reads.fa")):
    #     if line[0] != '>':
    #         reads.append(line.rstrip())
    # test for known genes
    aligner = Aligner(genome_sequence, genes)
    # alignments = []
    # hit_count = 0
    # start_time = time.time()
    # for r in reads:
    #     a = aligner.align(r)
    #     alignments.append(a)
    #     if a:
    #         hit_count += 1
    # print('read count: ' + str(len(reads)))
    # print('total time: ' + str(time.time() - start_time))
    # print('hit rate: ' + str(hit_count/len(reads)))

    aligner.align(read)

    # test findSubsequence
    # window = {((10, 11), (4, 5)), ((4, 7), (5, 8)), ((8, 9), (3, 4)), ((7, 8), (4, 5)), ((1, 2), (4, 5)), ((2, 4), (0, 2)), ((8, 9), (2, 3))}
    # sub = aligner.findRuns(window)
    # for seq in sub:
    #     print(seq)

def runKnownAndUnknown():
    with open('./genome.fa') as f:
        genome_sequence = f.readline()
        genome_sequence = f.readline() + '$'

    genes = set()
    isoforms = []
    exons = []
    for line in reversed(list(open("./genes.tab"))):
        elements = line.split('\t')
        if elements[0] == 'exon':
            ex = Exon(elements[1], int(elements[2]), int(elements[3]))
            exons.append(ex)
        elif elements[0] == 'isoform':
            iso = Isoform(elements[1], exons)
            isoforms.append(iso)
            exons = []
        elif elements[0] == 'gene':
            g = Gene(elements[1], isoforms)
            genes.add(g)
            isoforms = []

    reads = []
    for line in list(open("./reads.fa")):
        if line[0] != '>':
            reads.append(line.rstrip())

    aligner = Aligner(genome_sequence, genes)
    alignments = []
    hit_count = 0
    start_time = time.time()
    for r in reads:
        a = aligner.align(r)
        alignments.append(a)
        if a:
            hit_count += 1
    print('read count: ' + str(len(reads)))
    print('total time: ' + str(time.time() - start_time))
    print('known genes hit rate: ' + str(hit_count / len(reads)))

    genes = set()
    isoforms = []
    exons = []

    for line in reversed(list(open("./genes.tab"))):
        elements = line.split('\t')
        if elements[0] == 'unknown_exon':
            ex = Exon(elements[1], int(elements[2]), int(elements[3]))
            exons.append(ex)
        elif elements[0] == 'unknown_isoform':
            iso = Isoform(elements[1], exons)
            isoforms.append(iso)
            exons = []
        elif elements[0] == 'unknown_gene':
            g = Gene(elements[1], isoforms)
            genes.add(g)
            isoforms = []

    aligner = Aligner(genome_sequence, genes)
    alignments = []
    hit_count = 0
    start_time = time.time()
    for r in reads:
        a = aligner.align(r)
        alignments.append(a)
        if a:
            hit_count += 1
    print('read count: ' + str(len(reads)))
    print('total time: ' + str(time.time() - start_time))
    print('unknown genes hit rate: ' + str(hit_count / len(reads)))

def runFullVersion():
    with open('./genome.fa') as f:
        genome_sequence = f.readline()
        genome_sequence = f.readline() + '$'

    genes = set()
    isoforms = []
    exons = []
    for line in reversed(list(open("./genes.tab"))):
        elements = line.split('\t')
        if elements[0] == 'exon':
            ex = Exon(elements[1], int(elements[2]), int(elements[3]))
            exons.append(ex)
        elif elements[0] == 'isoform':
            iso = Isoform(elements[1], exons)
            isoforms.append(iso)
            exons = []
        elif elements[0] == 'gene':
            g = Gene(elements[1], isoforms)
            genes.add(g)
            isoforms = []

    reads = []
    for line in list(open("./reads.fa")):
        if line[0] != '>':
            reads.append(line.rstrip())

    aligner = Aligner(genome_sequence, genes)
    print('aligner initialized')
    alignments = []
    hit_count = 0
    start_time = time.time()
    subset_reads = 100

    for i in tqdm(range(subset_reads)):
        a = aligner.align(reads[i])
        alignments.append(a)
        if a:
            hit_count += 1
    print('read count: ' + str(len(alignments)))
    print('total time: ' + str(time.time() - start_time))
    print('total hit rate: ' + str(hit_count / len(alignments)))

def runSingleRead():
    read = 'ATTAGCGTGGCTCCTGTGCAAGGAAGACACACAAATTTGTGATGTGTTCC'
    with open('./genome.fa') as f:
        genome_sequence = f.readline()
        genome_sequence = f.readline() + '$'

    genes = set()
    isoforms = []
    exons = []
    for line in reversed(list(open("./genes.tab"))):
        elements = line.split('\t')
        if elements[0] == 'exon':
            ex = Exon(elements[1], int(elements[2]), int(elements[3]))
            exons.append(ex)
        elif elements[0] == 'isoform':
            iso = Isoform(elements[1], exons)
            isoforms.append(iso)
            exons = []
        elif elements[0] == 'gene':
            g = Gene(elements[1], isoforms)
            genes.add(g)
            isoforms = []

    aligner = Aligner(genome_sequence, genes)

    start_time = time.time()
    alignment = aligner.alignGenome(read)
    print(alignment)
    print('align time: ' + str(time.time() - start_time))
    print(read)

def testRadixSort():
    # s = 'ACGTAGCCG' * 2000 + '$'
    # s = 'ACTGGTTACCCTACTGATTAGGACTC$'
    # s = STRING
    # print(get_suffix_array(s) == naive_suffix_array(s))
    s = ''
    with open('./genome.fa') as f:
        s = f.readline()
        s = f.readline() + '$'
    get_suffix_array(s)

def runEval():
    known_genes = set()
    isoforms = []
    exons = []
    for line in reversed(list(open("./genes.tab"))):
        elements = line.split('\t')
        if elements[0] == 'exon':
            ex = Exon(elements[1], int(elements[2]), int(elements[3]))
            exons.append(ex)
        elif elements[0] == 'isoform':
            iso = Isoform(elements[1], exons)
            isoforms.append(iso)
            exons = []
        elif elements[0] == 'gene':
            g = Gene(elements[1], isoforms)
            known_genes.add(g)
            isoforms = []

    unknown_genes = set()
    isoforms = []
    exons = []
    for line in reversed(list(open("./genes.tab"))):
        elements = line.split('\t')
        if elements[0] == 'unknown_exon':
            ex = Exon(elements[1], int(elements[2]), int(elements[3]))
            exons.append(ex)
        elif elements[0] == 'unknown_isoform':
            iso = Isoform(elements[1], exons)
            isoforms.append(iso)
            exons = []
        elif elements[0] == 'unknown_gene':
            g = Gene(elements[1], isoforms)
            unknown_genes.add(g)
            isoforms = []

    all_known_isoforms = set()
    for k in known_genes:
        for i in k.isoforms:
            all_known_isoforms.add(i)

    all_unknown_isoforms = set()
    for u in unknown_genes:
        for i in u.isoforms:
            all_unknown_isoforms.add(i)

    with open('./genome.fa') as f:
        genome_sequence = f.readline()
        genome_sequence = f.readline() + '$'

    reads = []
    for line in list(open("./reads.fa")):
        if line[0] != '>':
            reads.append(line.rstrip())

    aligner = Aligner(genome_sequence, known_genes)
    print('aligner initialized')
    alignments = []
    hit_count = 0
    start_time = time.time()
    subset_reads = 100

    known_count = 0
    hidden_count = 0
    unaligned_count = 0
    index = index_isoform_locations(all_known_isoforms, all_unknown_isoforms)
    for i in range(subset_reads):
        read = reads[i]
        a = aligner.align(read)
        print(a)
        for m in a:
            read_seg = read[m[0]:m[0] + m[2]]
            ref_seg = genome_sequence[m[1]:m[1] + m[2]]
            # print(read_seg)
            # print(ref_seg)
        eval = evaluate_alignment(genome_sequence, reads[i], a, all_unknown_isoforms, index)
        if eval[0] == 'gene':
            known_count += 1
        if eval[0] == 'unaligned':
            unaligned_count += 1
        if eval[0] == 'hidden':
            hidden_count += 1
        alignments.append(a)
        if a:
            hit_count += 1
    print('known: ' + str(known_count))
    print('hidden: ' + str(hidden_count))
    print('unaligned: ' + str(unaligned_count))
    print('read count: ' + str(len(alignments)))
    print('total time: ' + str(time.time() - start_time))
    print('total hit rate: ' + str(hit_count / len(alignments)))

# testRadixSort()
# testAlignerInit()
# testBWTFunctions()
# runKnownAndUnknown()
# runFullVersion()
# runSingleRead()
runEval()
