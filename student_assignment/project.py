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
import math
import time
from shared import *

ALPHABET = [TERMINATOR] + BASES
MIN_INTRON_SIZE = 20
MAX_INTRON_SIZE = 10000
ANCHOR_LIMIT = 20

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
        # BWT-related structures for the reverse genome
        self.genome_sequence = genome_sequence
        self.reverse_genome = self.genome_sequence[::-1] + '$'
        self.reverse_sa = get_suffix_array(self.reverse_genome)
        self.reverse_bwt = get_bwt(self.reverse_genome, self.reverse_sa)
        self.reverse_F = get_F(self.reverse_bwt)
        self.reverse_M = get_M(self.reverse_F)
        self.reverse_occ = get_occ(self.reverse_bwt)
        
        self.known_genes = known_genes
        
        # BWT-related structures for reverse isoforms
        self.exon_indices = {}
        self.isoform_sequences = {}
        self.isoform_lengths = {}
        self.isoform_index_map = {}
        self.isoform_sas = {}
        self.isoform_M = {}
        self.isoform_occ = {}
        self.processIsoforms()
    
    def processIsoforms(self):
        """
        Populates 
            self.exon_indices
            self.isoform_lengths
            self.isoform_index_map
            self.isoform_sas
            self.isoform_M
            self.isoform_occ
        """
        for gene in self.known_genes:
            for isoform in gene.isoforms:
                i_id = isoform.id
                # the map_indices array is needed to keep track of where in the genome
                # a given MMP maps to, since self.mms() will return matches to the
                # concatenated sequence passed into the method
                map_indices = []
                self.exon_indices[i_id] = []
                sequence = ''
                for exon in sorted(isoform.exons, key=lambda item:item.start):
                    start = exon.start
                    end = exon.end
                    self.exon_indices[i_id].append((start, end))
                    map_indices.extend([i for i in range(start, end)])
                    sequence += self.genome_sequence[start:end]
                self.isoform_sequences[i_id] = sequence
                self.isoform_lengths[i_id] = len(sequence)
                self.isoform_index_map[i_id] = map_indices
                reverse_sequence = sequence[::-1] + '$'
                self.isoform_sas[i_id] = get_suffix_array(reverse_sequence)
                bwt = get_bwt(reverse_sequence, self.isoform_sas[i_id])
                F = get_F(bwt)
                self.isoform_M[i_id] = get_M(F)
                self.isoform_occ[i_id] = get_occ(bwt)   
        return

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
        
    def alignKnown(self, read_sequence):
        """
        Returns an alignment to known isoforms. 
        """
        # reverse read for finding MMPs using exact suffix match
        reverse_read = read_sequence[::-1]
        read_length = len(read_sequence)
        
        # seeds/anchor_seeds are a dictionaries mapping isoform ids to lists of
        # seeds and anchor seeds identified in the isoform (indexing isoform, not genome)
        # isoform_seeds = {}
        # isoform_windows = {}
        best_alignment = (None, None) # tuple where first value is isoform id and second value is alignment
        best_score = -math.inf
        for isoform in self.exon_indices.keys():
            isoform_sequence = self.isoform_sequences[isoform]
            seeds, anchor_seeds = self.findSeeds(read_sequence, isoform)
            windows = self.findWindows(seeds, anchor_seeds)
            for w in windows:
                runs = self.findRuns(w)
                for a in runs:
                    print(str(a) + ' len: ' + str(len(a)))
                    if a[0][1][0] != 0 or a[-1][1][1] != len(read_sequence):
                        continue
                    score, alignment = self.findAlignment(read_sequence, a, isoform_sequence)
                    if score > best_score:
                        best_score = score
                        best_alignment = (isoform, alignment)
        # print(self.formatAlignment(best_alignment[1], best_alignment[0])) 

    def findSeeds(self, read_sequence, isoform_id):
        """
        Finds seeds and marked anchor seeds for a given isoform id.
        """
        reverse_read = read_sequence[::-1]
        read_length = len(read_sequence)
        n = self.isoform_lengths[isoform_id]
        
        reverse_seeds = self.mms(reverse_read, read_length, self.isoform_M[isoform_id], self.isoform_occ[isoform_id])
        seeds = []
        anchor_seeds = [] # all the alignments that map less than 20 times are selected as anchors
        for match in reverse_seeds:
            genome_match, read_match = match
            g_start, g_end = genome_match
            r_start, r_end = read_match
            
            for i in range(g_start, g_end):
                reverse_start = self.isoform_sas[isoform_id][i]
                
                o_start = n - (reverse_start + (r_end - r_start))
                o_end = o_start + (r_end - r_start)
                genome_start = self.isoform_index_map[isoform_id][o_start]
                if o_end >= len(self.isoform_index_map[isoform_id]):
                    genome_end = self.isoform_index_map[isoform_id][o_end - 1] + 1
                else:
                    genome_end = self.isoform_index_map[isoform_id][o_end]
                
                original_r_end = read_length - r_start
                original_r_start = read_length - r_end
                
                # if there is an exon junction here, create separate seeds
                if genome_end - genome_start > o_end - o_start:
                    last_end = o_start
                    last_r_end = original_r_start
                    count = 0
                    for i in range(o_start, o_end - 1):
                        count += 1
                        if self.isoform_index_map[isoform_id][i] < self.isoform_index_map[isoform_id][i + 1] - 1:
                            seeds.append(((last_end, i + 1), (last_r_end, last_r_end + count)))
                            if g_end - g_start < ANCHOR_LIMIT:
                                seeds.append(((last_end, i + 1), (last_r_end, last_r_end + count)))
                            last_r_end = last_r_end + count
                            last_end = i + 1
                            count = 0
                    seeds.append(((last_end, o_end), (last_r_end, original_r_end)))
                    if g_end - g_start < ANCHOR_LIMIT:
                        seeds.append(((last_end, o_end), (last_r_end, original_r_end)))
                else:
                    seeds.append(((o_start, o_end), (original_r_start, original_r_end)))
                    if g_end - g_start < ANCHOR_LIMIT:
                        anchor_seeds.append(((o_start, o_end), (original_r_start, original_r_end)))
        return seeds, anchor_seeds

    def alignGenome(self, read_sequence):
        """
        Returns an alignment to known isoforms. Separate from alignIsoforms since we
        constructed BWT-related structures in the __init__ and don't want to repeat those operations
        """
        reverse_seeds = self.mms(read_sequence[::-1], len(read_sequence), self.reverse_M, self.reverse_occ)
        seeds = []
        anchor_seeds = []
        n = len(self.genome_sequence)
        read_length = len(read_sequence)
        for match in reverse_seeds:
            genome_match, read_match = match
            g_start, g_end = genome_match
            r_start, r_end = read_match
            for i in range(g_start, g_end):
                reverse_start = self.reverse_sa[i]
                o_start = n - (reverse_start + (r_end - r_start))
                o_end = o_start + (r_end - r_start)
                original_r_end = read_length - r_start
                original_r_start = read_length - r_end
                seed = ((o_start, o_end), (original_r_start, original_r_end))
                seeds.append(seed)
                if o_start - o_end < ANCHOR_LIMIT:
                    anchor_seeds.append(seed)
        best_alignment = None
        best_score = -math.inf
        windows = self.findWindows(seeds, anchor_seeds)
        for w in windows:
            runs = self.findRuns(w)
            for a in runs:
                score, alignment = self.findAlignment(read_sequence, a)
                if score > best_score:
                    best_score = score
                    best_alignment = alignment
        return self.formatAlignment(best_alignment)         
        
    def mms(self, read, i, M, occ):
        """
        Returns a set of tuples where each element is a ((sa_start, sa_end), (start, end)) 
        tuple representing start and end indices of the read in the 
        genome, and the chunk of the read we're using
        
        Input:
            read: the pattern string
            i: consider read[:i] for exact_suffix_matches
            M: M array for use by exact_suffix_matches
            occ: occ array for use by exact_suffix_matches
        """
        if i <= 0:
            return {}
        max_suffix = exact_suffix_matches(read[:i], M, occ)
        sa_indices, length = max_suffix
        if sa_indices == None:
            return {}
        start = sa_indices[0]
        end = sa_indices[1]
        return {((start, end), (i - length, i))}.union(self.mms(read, i - length, M, occ))

    def findWindows(self, seeds, anchor_seeds, window_size=MAX_INTRON_SIZE):
        """
            Return genomic windows in a 2D array of seeds
        """
        windows = []
        for a in anchor_seeds:
            window = []
            # seeds are in the format ((i_s, i_e), (r_s, r_e))
            window_start = a[0][0] - window_size
            if window_start < 0:
                window_start = 0
            window_end = a[0][1] + window_size
            for s in seeds:
                if s[0][0] >= window_start and s[0][1] <= window_end:
                    window.append(s)
            if window not in windows:
                windows.append(window)
        return windows
    
    def findRuns(self, window):
        """
            Return list of sets of seeds in the windows, where seeds
            are strictly increasing. We return the longest sets.
            
            Increasing = order is preserved between read and genome. If
            there are multiple sequences that are the same length, return all of them.
        """
        w = sorted(window, key=lambda item: item[0][0])
        store = []
        def increasingSubsequence(sorted_window, output):
            nonlocal store
            if len(sorted_window) == 0: 
                if len(output) != 0: 
                    # storing result 
                    store.append(output) 
                return
            increasingSubsequence(sorted_window[1:], output[:])
            if len(output) == 0: 
                increasingSubsequence(sorted_window[1:], sorted_window[:1]) 
            elif sorted_window[0][1][0] > output[-1][1][0] and sorted_window[0][0][0] > output[-1][0][0]: 
                output.append(sorted_window[0]) 
                increasingSubsequence(sorted_window[1:], output[:]) 
        increasingSubsequence(w, [])
        max_len = len(max(store, key=len))
        seen = set()
        return [s for s in store if len(s) == max_len and tuple(s) not in seen and not seen.add(tuple(s))]
    
    def findAlignment(self, read_sequence, seeds, isoform_sequence=None):
        """
            Returns best alignment and score
        """
        best_alignment = []
        total_score = 0
        reference_sequence = None
        # If no sequence passed in, assume alignment to genome
        if isoform_sequence == None:
            genome_start = seeds[0][0][0]
            genome_end = seeds[-1][0][1]
            reference_sequence = self.genome_sequence[genome_start, genome_end]
        else:
            reference_sequence = isoform_sequence
        # for every intron in the read
        for i in range(len(seeds) - 1):
            best_alignment.append(seeds[i]) # all seeds make it to the best alignment
            total_score += seeds[i][1][1] - seeds[i][1][0]
            best_indel_index = 0
            best_indel_score = 0
            read_intron = read_sequence[seeds[i][1][1]:seeds[i + 1][1][0]] # end of first seed to start of next
            read_intron_length = len(read_intron)
            reference_intron = reference_sequence[seeds[i][0][1]:seeds[i + 1][0][0]] # end of first seed to start of next
            reference_intron_length = len(reference_intron)
            if read_intron_length == 0: # only consider gap penalty
                total_score += (seeds[i + 1][0][0] - seeds[i][0][1]) * -1
                continue
            else:
                if read_intron_length < reference_intron_length:
                    for j in range(read_intron_length): # for every indel location
                        first_read_seg = read_intron[:j]
                        last_read_seg = read_intron[j:]
                        first_iso_seg = reference_intron[:j]
                        last_iso_seg = reference_intron[reference_intron_length - (read_intron_length - j):]
                        gap_penalty = (reference_intron_length - read_intron_length) * -1
                        score = self.matchScore(first_read_seg, first_iso_seg) + self.matchScore(last_read_seg, last_iso_seg) + gap_penalty
                        if score > best_indel_score:
                            best_indel_index = j
                            best_indel_score = 0
                    read_part_1 = read_sequence[seeds[i][1][1]:seeds[i][1][1] + best_indel_index]
                    read_part_2 = read_sequence[seeds[i][1][1] + best_indel_index:seeds[i + 1][1][0]]
                    seed1 = ((seeds[i][0][1], seeds[i][0][1] + len(read_part_1)), (seeds[i][1][1], seeds[i][1][1] + best_indel_index))
                    seed2 = ((seeds[i + 1][0][0] - len(read_part_2), seeds[i + 1][0][0]), (seeds[i][1][1] + best_indel_index, seeds[i + 1][1][0]))
                    best_alignment.append(seed1)
                    best_alignment.append(seed2)
                    total_score += best_indel_score
                else:
                    # If the read intron is longer than the reference intron, this is not a valid alignment (introduces gaps)
                    return 0, []
        return total_score, best_alignment
    
    def matchScore(self, seq1, seq2):
        """
            Scores the alignment between seq1 and seq2, assumed to have the same length
        """
        if len(seq1) != len(seq2):
            raise ValueError('matchScore only takes string inputs of equal length')
        score = 0
        for i in range(len(seq1)):
            if seq1[i] == seq2[i]:
                score += 1
            else:
                score -= 1
        return score

    def calculateSeedsLength(self, seeds):
        """
            Helper function to return length of all seeds 
        """
        length = 0
        for s in seeds:
            length += s[1][1] - s[1][0]
        return length
            
    
    def formatAlignment(self, seeds, isoform_id=None):
        """
            Format alignment for project output. If isoform_id is passed in, we need to convert
            indices back to genome indices. If not, we keep these indices.
            
            "We expect you to specify your alignment as a python list of k tuples of (start index, genome start
            index, length)"
        """
        alignment = []
        if seeds == None:
            return []
        if isoform_id == None:
            for s in seeds:
                alignment.append((s[1][0], s[0][0], s[1][1] - s[1][0]))
        else:
            for s in seeds:
                genome_start = self.isoform_index_map[isoform_id][s[0][0]]
                alignment.append((s[1][0], genome_start, s[1][1] - s[1][0]))
        return alignment
        