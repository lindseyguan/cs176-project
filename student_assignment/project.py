""" 
    RNA Alignment Assignment
    
    Implement each of the functions below using the algorithms covered in class.
    You can construct additional functions and data structures but you should not
    change the functions' APIs.

    You will be graded on the helper function implementations as well as the RNA alignment, although
    you do not have to use your helper function.
    
    *** Make sure to comment out any print statement so as not to interfere with the grading script
"""

import sys # DO NOT EDIT THIS
import time # DELETE THIS LATER
import numpy as np
from shared import *

ALPHABET = [TERMINATOR] + BASES
radix_k = 4

class Bucket:
    def __init__(self, bucket_id, str_arr, k):
        self.bucket_id = bucket_id
        self.str_arr = str_arr
        self.k = k
        self.sub_buckets = {}
    
    def get_sub_buckets(self):
        """
        Recursively generates sub buckets using the kth character of the strings
        """
        if self.k >= radix_k:
            return {}
        
        buckets = {}
        for p in self.str_arr:
            char = p[1][self.k]
            if char in buckets:
                buckets[char].str_arr.append(p)
            else:
                new_bucket = Bucket(char, [p], self.k + 1)
                buckets[char] = new_bucket
        for b in buckets.values():
            b.get_sub_buckets()
            
        self.sub_buckets = buckets
    
def lex_traverse(bucket, s):
    """
    Recursively returns a tuple of all indices of all strings of all sub buckets of a bucket 
    in lexicographic order
    """
    traversed = []
    if bucket.str_arr == []:
        return []
    if bucket.sub_buckets == {}:
        return sorted(bucket.str_arr, key=lambda item: s[item[0]:])
    else:
        for key in sorted(bucket.sub_buckets):
            traversed.extend(lex_traverse(bucket.sub_buckets[key], s))
    return traversed
        
        
def naive_suffix_array(s):
    start_time = time.time()
    index_suffix_dict = {i:s[i:] for i in range(len(s))}
    a = [k for k, v in sorted(index_suffix_dict.items(), key=lambda item: item[1])]
    print('naive: ' + str((time.time() - start_time) * 1000))
    # print(a)
    return a
    
def get_suffix_array(s):
    """
    Naive implementation of suffix array generation (0-indexed). You do not have to implement the
    KS Algorithm. Make this code fast enough so you have enough time in Aligner.__init__ (see bottom).

    Input:
        s: a string of the alphabet ['A', 'C', 'G', 'T'] already terminated by a unique delimiter '$'
    
    Output: list of indices representing the suffix array

    >>> get_suffix_array('GATAGACA$')
    [8, 7, 5, 3, 1, 6, 4, 0, 2]
    """
    start_time = time.time()
    str_arr = []
    n = len(s)
    for i in range(n):
        if i + radix_k + 1> n:
            str_padded = s[i:].ljust(radix_k + 1)
        else:
            str_padded = s[i:i+radix_k+1]
        str_arr.append((i, str_padded))
    main_bucket = Bucket('MAIN', str_arr, 0)
    main_bucket.get_sub_buckets()
    radix_sorted = [int(x[0]) for x in lex_traverse(main_bucket, s)]
    print('radix: ' + str((time.time() - start_time) * 1000))
    # print(radix_sorted)
    # print(naive_suffix_array(s))
    print(radix_sorted == naive_suffix_array(s))


def get_bwt(s, sa):
    """
    Input:
        s: a string terminated by a unique delimiter '$'
        sa: the suffix array of s

    Output:
        L: BWT of s as a string
    """
    L = ''
    n = len(sa)
    for i in sa:
        L += s[(i + n - 1) % n]
    return L

def get_F(L):
    """
    Input: L = get_bwt(s)

    Output: F, first column in Pi_sorted
    """
    return ''.join(sorted(L))

def get_M(F):
    """
    Returns the helper data structure M (using the notation from class). M is a dictionary that maps character
    strings to start indices. i.e. M[c] is the first occurrence of "c" in F.

    If a character "c" does not exist in F, you may set M[c] = -1
    """
    M = {}
    for c in ALPHABET:
        if c not in M:
            M[c] = F.find(c)
    return M

def get_occ(L):
    """
    Returns the helper data structure OCC (using the notation from class). OCC should be a dictionary that maps 
    string character to a list of integers. If c is a string character and i is an integer, then OCC[c][i] gives
    the number of occurrences of character "c" in the bwt string up to and including index i
    """
    occ = {}
    for c in ALPHABET:
        occ[c] = []
    for i in range(len(L)):
        for k in occ.keys():
            if L[i] == k:
                if i == 0:
                    occ[k] += [1]
                else:
                    occ[k] += [occ[k][i - 1] + 1]
            else:
                if i == 0:
                    occ[k] += [0]
                else:
                    occ[k] += [occ[k][i - 1]]
    return occ

def exact_suffix_matches(p, M, occ):
    """
    Find the positions within the suffix array sa of the longest possible suffix of p 
    that is a substring of s (the original string).
    
    Note that such positions must be consecutive, so we want the range of positions.

    Input:
        p: the pattern string
        M, occ: buckets and repeats information used by sp, ep

    Output: a tuple (range, length)
        range: a tuple (start inclusive, end exclusive) of the indices in sa that contains
            the longest suffix of p as a prefix. range=None if no indices matches any suffix of p
        length: length of the longest suffix of p found in s. length=0 if no indices matches any suffix of p

        An example return value would be ((2, 5), 7). This means that p[len(p) - 7 : len(p)] is
        found in s and matches positions 2, 3, and 4 in the suffix array.

    >>> s = 'ACGT' * 10 + '$'
    >>> sa = get_suffix_array(s)
    >>> sa
    [40, 36, 32, 28, 24, 20, 16, 12, 8, 4, 0, 37, 33, 29, 25, 21, 17, 13, 9, 5, 1, 38, 34, 30, 26, 22, 18, 14, 10, 6, 2, 39, 35, 31, 27, 23, 19, 15, 11, 7, 3]
    >>> L = get_bwt(s, sa)
    >>> L
    'TTTTTTTTTT$AAAAAAAAAACCCCCCCCCCGGGGGGGGGG'
    >>> F = get_F(L)
    >>> F
    '$AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT'
    >>> M = get_M(F)
    >>> sorted(M.items())
    [('$', 0), ('A', 1), ('C', 11), ('G', 21), ('T', 31)]
    >>> occ = get_occ(L)
    >>> type(occ) == dict, type(occ['$']) == list, type(occ['$'][0]) == int
    (True, True, True)
    >>> occ['$']
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    >>> exact_suffix_matches('ACTGA', M, occ)
    ((1, 11), 1)
    >>> exact_suffix_matches('$', M, occ)
    ((0, 1), 1)
    >>> exact_suffix_matches('AA', M, occ)
    ((1, 11), 1)
    """
    # COUNT algorithm as described in Lec 7
    c = p[-1]
    sp = M[c]
    M_sorted = [[k, v] for k, v in sorted(M.items(), key=lambda item: item[1])]
    if not M_sorted.index([c, M[c]]) >= len(M_sorted) - 1:
        ep = M[M_sorted[M_sorted.index([c, M[c]]) + 1][0]] - 1
    else:
        ep = len(occ[ALPHABET[0]]) - 1
    i = len(p) - 1
    while i > 0:
        i -= 1
        old_sp = sp
        old_ep = ep
        sp = M[p[i]] + occ[p[i]][sp - 1]
        ep = M[p[i]] + occ[p[i]][ep] - 1
        if ep <= sp:
            i += 1
            sp = old_sp
            ep = old_ep
            break
    return ((sp, ep + 1), len(p) - i)


MIN_INTRON_SIZE = 20
MAX_INTRON_SIZE = 10000

class Aligner:
    def __init__(self, genome_sequence, known_genes):
        """
        Initializes the aligner. Do all time intensive set up here. i.e. build suffix array.

        genome_sequence: a string (NOT TERMINATED BY '$') representing the bases of the of the genome
        known_genes: a python set of Gene objects (see shared.py) that represent known genes. You can get the isoforms 
                     and exons from a Gene object

        Time limit: 500 seconds maximum on the provided data. Note that our server is probably faster than your machine, 
                    so don't stress if you are close. Server is 1.25 times faster than the i7 CPU on my computer

        """
        self.sa = get_suffix_array(genome_sequence)
        self.known_genes = known_genes

    def align(self, read_sequence):
        """
        Returns an alignment to the genome sequence. An alignment is a list of pieces. 
        Each piece consists of a start index in the read, a start index in the genome, and a length 
        indicating how many bases are aligned in this piece. Note that mismatches are count as "aligned".

        Note that <read_start_2> >= <read_start_1> + <length_1>. If your algorithm produces an alignment that 
        violates this, we will remove pieces from your alignment arbitrarily until consecutive pieces 
        satisfy <read_start_2> >= <read_start_1> + <length_1>

        Return value must be in the form (also see the project pdf):
        [(<read_start_1>, <reference_start_1, length_1), (<read_start_2>, <reference_start_2, length_2), ...]

        If no good matches are found: return the best match you can find or return []

        Time limit: 0.5 seconds per read on average on the provided data.
        """
        pass

def testBWTFunctions():
    # s = 'ACGT' * 10 + '$'
    # # print(s)
    # sa = get_suffix_array(s)
    # # print(sa)
    # L = get_bwt(s, sa)
    # # print(L)
    # F = get_F(L)
    # # print(F)
    # M = get_M(F)
    # # print(M)
    # occ = get_occ(L)
    # # print(occ)
    # matches = exact_suffix_matches('$', M, occ)
    # # print(matches)
    pass

def testAlignerInit():
    # Testing runtime of Aligner init
    genome_sequence = ''
    with open('./genome.fa') as f:
        genome_sequence = f.readline()
        genome_sequence = f.readline() + '$'

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

    start_time = time.time()
    aligner = Aligner(genome_sequence, genes)
    print(time.time() - start_time)

def testRadixSort():
    s = 'ACGTAGCCG' * 20000 + '$'
    # s = 'ACGACGACG$'
    naive_suffix_array(s)
    get_suffix_array(s)

testRadixSort()
