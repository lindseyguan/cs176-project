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
import ctypes
import ctypes.util
import functools
import numpy as np
import time
from shared import *

ALPHABET = [TERMINATOR] + BASES

libc_name = ctypes.util.find_library("c")
libc = ctypes.CDLL(libc_name)
libc.memcmp.argtypes = (
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.c_size_t,
)

def sort_substrings(suffixes, string_ptr):
    def cmp_func(a, b):
        a_start, a_end = a
        b_start, b_end = b
        a_len = a_end - a_start
        b_len = b_end - b_start
        min_len = a_len
        if b_len < a_len:
            min_len = b_len

        # call C's memcmp (!)
        cmp_value = libc.memcmp(string_ptr + a_start, string_ptr + b_start, min_len)

        if cmp_value != 0:
            return cmp_value
        else:
            # whichever string is shorter
            return b_len - a_len

    suffixes.sort(key=functools.cmp_to_key(cmp_func))
    return suffixes
    
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
    suffixes = [(i, len(s)) for i in range(0, len(s))]
    s = s.encode()
    string_bb = bytearray(s)
    string_p = ctypes.cast(s, ctypes.c_char_p)
    string_p2 = (ctypes.c_char * len(s)).from_buffer(string_bb)
    string_ptr = ctypes.cast(string_p2, ctypes.c_void_p).value
    sorted_suffixes = sort_substrings(suffixes, string_ptr)
    # print("memcmp sort: " + str((time.time() - start_time) * 1000))
    return [s[0] for s in sorted_suffixes]

def naive_suffix_array(s):
    """
    Naive implementation of suffix array generation, only used to check correctness
    """
    start_time = time.time()
    index_suffix_dict = {i:s[i:] for i in range(len(s))}
    a = [k for k, v in sorted(index_suffix_dict.items(), key=lambda item: item[1])]
    print('naive: ' + str((time.time() - start_time) * 1000))
    # print(a)
    return a

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
        if ep < sp:
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
        self.genome_sequence = '$' + genome_sequence
        print('genome: ' + self.genome_sequence)
        self.reverse_genome = self.genome_sequence[::-1]
        print('reverse genome: ' + self.reverse_genome)
        self.reverse_sa = get_suffix_array(self.reverse_genome)
        self.reverse_bwt = get_bwt(self.reverse_genome, self.reverse_sa)
        self.reverse_F = get_F(self.reverse_bwt)
        self.reverse_M = get_M(self.reverse_F)
        self.reverse_occ = get_occ(self.reverse_bwt)
        
        self.known_genes = known_genes
        self.isoform_indices = self.processIsoforms()

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
    
    def processIsoforms(self):
        """
        Returns dictionary of known gene indices grouped by isoform
        
        Output:
            indices: dictionary of lists, where each element is (start, end) in the genome
        """
        genes = {}
        for gene in self.known_genes:
            for isoform in gene.isoforms:
                id = isoform.id
                genes[id] = {}
                for exon in sorted(isoform.exons, key=lambda item:item.start):
                    genes[id] = (exon.start, exon.end)
        return genes
                
    def mms(self, read, i, mapKnownGenes=True):
        """
        Returns a set of tuples where each element is a ((sa_start, sa_end), (start, end)) 
        tuple representing start and end indices of the read in the 
        genome, and the chunk of the read we're using
        
        Input:
            read: the pattern string
            i: consider read[:i] for exact_suffix_matches
            M: M array to be used for suffix matching
            occ: occ array to be used for suffix matching
            mapKnownGenes: if true, map to known genes. if false, map to genome
        """
        if i <= 0:
            return {}
        pattern = read[:i]
        max_suffix = exact_suffix_matches(pattern, self.reverse_M, self.reverse_occ)
        sa_indices, length = max_suffix
        if sa_indices == None:
            return {}
        start = sa_indices[0]
        end = sa_indices[1]
        # FIX ME
        return {((start, end), (i - length, i))}.union(self.mmp(read, i - length, mapKnownGenes))

