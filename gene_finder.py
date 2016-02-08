# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Carl Moser

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
import time


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    s = ''
    for c in dna:
        s = get_complement(c) + s
    return s


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGATGG")
    'ATGAGATGG'
    """
    stops = {'TAA', 'TAG', 'TGA'}
    i = 0
    while i <= len(dna):
        cod = dna[i:i+3]
        if cod in stops:
            return dna[0:i]
        i += 3
    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGAATGTGCCC")
    ['ATGCATGAATGTAGAATGTGCCC']
    """
    orfs = []
    index = len(dna)//3
    prevorf = ''
    for i in range(index):
        if(dna[3*i:3*i+3] == 'ATG'):
            toappend = rest_of_ORF(dna[3*i:])
            if toappend not in prevorf:
                orfs.append(toappend)
                prevorf = toappend
    return orfs

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    frames = []
    for i in range(3):
        for orf in find_all_ORFs_oneframe(dna):
            frames.append(orf)
        dna = dna[1:]
    return frames



def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    >>> find_all_ORFs_both_strands("AAAA")
    []
    """
    orfs = find_all_ORFs(dna)
    for orf in find_all_ORFs(get_reverse_complement(dna)):
        orfs.append(orf)
    return orfs

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("")
    ''
    """
    orfs = find_all_ORFs_both_strands(dna)
    if len(orfs) == 0:
        return ''
    longest = orfs[0]
    for orf in range(len(orfs)):
        if len(orfs[orf]) > len(longest):
            longest = orfs[orf]
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longest = 0
    for trial in range(num_trials):
        curval = len(longest_ORF(shuffle_string(dna)))
        if curval > longest:
            longest = curval
    return longest



def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    strand = ''
    i = 0
    while i < len(dna)//3:
        cod = dna[3*i:3*i+3]
        strand += aa_table[cod]
        i += 1
    return strand


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    print(threshold)
    strands = find_all_ORFs_both_strands(dna)
    longer = []
    for strand in strands:
        if len(strand) >= threshold:
            coded = coding_strand_to_AA(strand)
            longer.append(coded)
    return longer


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    dna = load_seq("./data/X73525.fa")
    start_time = time.time()
    genes = gene_finder(dna)
    print("--- %s seconds ---" % (time.time() - start_time))
    print(genes)