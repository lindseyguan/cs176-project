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
    s = 'ACTGGTTACCCTACTGATTAGGACTC$'
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
    matches = exact_suffix_matches('CTC', M, occ)
    print(matches)

def testAlignerInit():
    # Testing runtime of Aligner init
    # genome_sequence = ''
    # with open('./genome.fa') as f:
    #     genome_sequence = f.readline()
    #     genome_sequence = f.readline() + '$'
    genome_sequence = 'ACTGGTTACCCTACTGATTAGGACTC'
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
    start_time = time.time()
    aligner = Aligner(genome_sequence, genes)
    read = 'TACCG'
    n = len(genome_sequence)
    r_n = len(read)
    mms = aligner.mms(read[::-1], r_n)
    print(mms)
    sa_naive = naive_suffix_array(genome_sequence[::-1] + '$')
    # for i in range(n):
    #     print(str(i) + ': ' + str(sa_naive[i]) + ' -> ' + genome_sequence[sa_naive[i]:n])
    for match in mms:
        genome_match, read_match = match
        g_start, g_end = genome_match
        r_start, r_end = read_match
        print('pattern: ' + read[::-1][r_start:r_end])
        for i in range(g_start, g_end):
            print('match: ' + str(sa_naive[i]) + ', or ' + genome_sequence[::-1][sa_naive[i]:])
        
    print(time.time() - start_time)

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
