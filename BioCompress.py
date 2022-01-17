import numpy as np
from SuffixTree import *


"""Recursive function for calculating the representation
  of the number decnum in base radix."""
def d2r(radix, decnum, current):
    res = decnum % radix
    decnum = decnum // radix
    if decnum == 0:
        return str(res) + current
    else:
        return d2r(radix, decnum, str(res) + current)

"""Wrapper function for d2r."""
def dec2radix(radix, decnum):
    return d2r(radix, decnum, '')

"""Looks for the longest string that starts at the 
   search buffer and that is a prefix of the lookahead
   buffer, using a SUFFIX TREE for the search. This method
   traverse the SuffixTree exploting the propertie of
   McCreights' SuffixTree construction: only the first symbol
   of the edge label must be checked in every node."""
def ST_based_reproducible_extension(STD, s, lim=None, limsup=None):
    v = STD.root
    i = 0
    l = 0
    while i < len(s):
        if v.has_transition(s[i]):
            if lim is not None and v.transition_links[s[i]].new_idx < lim:
                return
            if limsup is not None and v.transition_links[s[i]].idx >= limsup:
                return (v.idx, l)
            label = STD.edgeLabel(v.transition_links[s[i]], v)
            idx = 0
            while i+idx < len(s) and idx < len(label):
                if s[i+idx] == label[idx]:
                    idx = idx + 1
                    l = l + 1
                else:
                    if lim is not None:
                        return (v.transition_links[s[i]].new_idx, l)
                    else:
                        return (v.transition_links[s[i]].idx, l)
            v = v.transition_links[s[i]]
            i = i + idx
        else:
            if lim is not None:
                return (v.new_idx, l)
            else:
                return (v.idx, l)
    if lim is not None:
        return (v.new_idx, l)
    else:
        return (v.idx, l)

"""Convert list to text string"""
def list_to_string(L):
    S = ""
    for l in L:
        S += str(l)
    return S

"""Get the palindrome: chaging each symbol
   with it complementary and reversing this string"""
def palindrome(S):
    SR = ""
    for idx, s in enumerate(S):
        if s == 'A':
            SR += 'T'
        elif s == 'T':
            SR += 'A'
        elif s == 'C':
            SR += 'G'
        elif s == 'G':
            SR += 'C'
    return SR[::-1]

"""radix representation of a decimal number in 'pad' bits.
   Padding is done if necessary."""
def encode_padbits(decnum, radix, pad):
    s = dec2radix(radix, decnum)
    #print("radix rep:= ", s)
    if len(s) < pad:
        return '0'*(pad - len(s)) + s
    return s

"""return n-th fibonacci number calculation"""
def fibonacci(n):
    a, b = 0, 1
    for i in range(0, n):
        a, b = b, a + b
    return a

"""Perform fibonacci codeword generation."""
def fibonacci_encode_number(n):
    code = [0]
    
    remind = n
    while remind != 0:
        idx = 1
        fib = fibonacci(idx)
        while fibonacci(idx) <= remind:
            fib = fibonacci(idx)
            idx += 1
            code.append(0)
        code[idx-2] = 1
        remind = remind-fib
    return code+[1]

def encode_pair(par):
    x = par[0]
    y = par[1]

    if x < 7:
        cod_x = encode_padbits(x, radix=2, pad=3)
    elif 7 <= x and x <= 38:
        cod_x = encode_padbits(x-7, radix=2, pad=5)
    else:
        cod_x = list_to_string(fibonacci_encode_number(x))

    if y < 7:
        cod_y = encode_padbits(y, radix=2, pad=3)
    elif 7 <= y and y <= 38:
        cod_y = encode_padbits(y-7, radix=2, pad=5)
    else:
        cod_y = list_to_string(fibonacci_encode_number(y))
    
    return (cod_x, cod_y)

def encode_symbol(simb):
    if simb == 'A':
        return '00'
    elif simb == 'C':
        return '01'
    elif simb == 'G':
        return '10'
    elif simb == 'T':
        return '11'

def BioCompress(D, S):
    totcode = []
    ST = STree(D)

    i = 0
    while i < len(S):
        pos, size = ST_based_reproducible_extension(ST, S[i:], limsup=i)
        pos_p, size_p = ST_based_reproducible_extension(ST, palindrome(S[i:]), limsup=i)

        if size <= 1:
            encoded_symbol = encode_symbol(S[i])
            #print("encoded symbol: ", encoded_symbol)
            totcode.append(encoded_symbol)
            size = 1
        else:
            encoded_pair = encode_pair((pos, size))
            #print("encoded pair: ", encoded_pair)
            totcode.append(encoded_pair)
        i = i+size

    print(totcode)

def BioCompress_with_palindromes(D, S, c):
    totcode = []
    ST = STree(D)

    i = 0
    while i < len(S):
        pos, size = ST_based_reproducible_extension(ST, S[i:], limsup=i)
        pos_p, size_p = ST_based_reproducible_extension(ST, palindrome(S[i:]), limsup=i)

        #print((pos, size), " vs ", (pos_p, size_p))
        encoded_pair = encode_pair((pos, size))
        encoded_pair_p = encode_pair((pos_p, size_p))
        #print("encoded pairs: ", encoded_pair, " vs ", encoded_pair_p)
        if len(encoded_pair[0]) + len(encoded_pair[1]) + c > 2*size:
            if len(encoded_pair_p[0]) + len(encoded_pair_p[1]) + c > 2*size_p:
                encoded_symbol = encode_symbol(S[i])
                #print("encoded symbol: ", encoded_symbol)
                totcode.append(encoded_symbol)
                i += 1
            else:
                #print("palindrome: ", encoded_pair_p)
                totcode.append(encoded_pair_p)
                i += size_p
        else:
                if len(encoded_pair[1]) > len(encoded_pair_p[1]) or len(encoded_pair_p[0]) + len(encoded_pair_p[1]) + c > 2*size_p:
                    #print("factor: ", encoded_pair)
                    totcode.append(encoded_pair)
                    i += size
                else:
                    #print("palindrome: ", encoded_pair_p)
                    totcode.append(encoded_pair_p)
                    i += size_p
    print(totcode)


cad = 'AACTGTTGTTGTTAGAACTGTTGTT'
print(cad)
print(len(cad))

BioCompress(D=cad, S=cad)
BioCompress_with_palindromes(D=cad, S=cad, c=0)
print("Almost DONE!")
