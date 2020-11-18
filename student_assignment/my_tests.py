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
from shared import *
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
    genome_sequence = ''
    with open('./genome.fa') as f:
        genome_sequence = f.readline()
        genome_sequence = f.readline() + '$'
    # genome_sequence = 'ACTGGTTACCCTACTGACCG'
    read = 'ATTACTCTTGGGAATGAAATCCTATCTATATAAGCTGTGGTTTGAAATCC'
    genes = set()

    gene_id = ''
    isoforms = []

    isoform_id = ''
    exons = []

    exon_id = ''
    start = 0
    end = 0
        
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

    # test MMP for known genes
    aligner = Aligner(genome_sequence, genes)
    print(aligner.alignGenome(read))
    # aligner.alignKnown(read)
    print(read)
    print(genome_sequence[10359306:10359306+14])
    
    # test findSubsequence
    # window = {((10, 11), (4, 5)), ((4, 7), (5, 8)), ((8, 9), (3, 4)), ((7, 8), (4, 5)), ((1, 2), (4, 5)), ((2, 4), (0, 2)), ((8, 9), (2, 3))}
    # sub = aligner.findRuns(window)
    # for seq in sub:
    #     print(seq)
    

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
    
# testRadixSort()
testAlignerInit()
# testBWTFunctions()
