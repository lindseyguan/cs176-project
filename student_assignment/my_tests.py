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
    # genome_sequence = ''
    # with open('./genome.fa') as f:
    #     genome_sequence = f.readline()
    #     genome_sequence = f.readline() + '$'
    genome_sequence = 'ACTGGTTACCCTACTGA'
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

    test_exons = [Exon('ENSE00001802701', 8250613, 8250877), Exon('ENSE00001729938', 8252369, 8252739)]
    test_isoforms = [Isoform('ENST00000433210', test_exons)]
    test_gene = Gene('ENSG00000231620', test_isoforms)
    assert(test_gene in genes)

    # test MMS
    aligner = Aligner(genome_sequence, genes)
    read = 'TACCTA'
    n = len(genome_sequence)
    r_n = len(read)
    aligner.align(read)

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
