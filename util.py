import re
import sys
import time

def pairs(s):
    '''Yields two things at a time, repeating "middle" things.'''
    while True:
        try:
            thing1 = s[0]
            thing2 = s[1]
            s      = s[1:]
        except IndexError:
            return
        yield thing1, thing2

####################################################

# Utility functions for taking reverse complements

def complement(string):
    '''
    Regex for substituting a base with its complement.
    A <--> T; C <--> G; N stays N.
    R (A or G) <--> Y (T or C); K (G or T) <--> M (C or A).
    W (A or T) stays W; S (G or C) stays S.
    H (not-G) <--> D (not-C); B (not-A) <--> V (not-T).
    http://www.bioinformatics.org/sms/iupac.html
    '''
    base = string.group(0).upper()
    if base == 'A': return 'T'
    elif base == 'T': return 'A'
    elif base == 'C': return 'G'
    elif base == 'G': return 'C'
    elif base == 'N': return 'N'
    elif base == 'R': return 'Y'
    elif base == 'Y': return 'R'
    elif base == 'S': return 'S'
    elif base == 'W': return 'W'
    elif base == 'K': return 'M'
    elif base == 'M': return 'K'
    elif base == 'H': return 'D'
    elif base == 'B': return 'V'
    elif base == 'V': return 'B'
    elif base == 'D': return 'H'
    else: print "unrecognized character"; print base; sys.exit(2)
    # TODO: better exit?

def reverse_complement(sequence):
    '''Returns the reverse complement of the input sequence.'''
    rev = sequence[::-1]        # reverses input string
    return re.sub('.', complement, rev)    # regex for substitute each base with its complement

####################################################

def printwithtime(arg):
    '''Wrapper around "print"'''
    print "[" + time.strftime('%X %x %Z') + "] " + str(arg)
