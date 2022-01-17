#Implementation of the LZ77 compression algorithm.

import numpy as np
from SuffixTree import *

#Coder----------------------------------------------------------------
class Coder:
    def __init__(self, alphabet, n, Ls):
        self.alph = alphabet
        self.n = n
        self.Ls = Ls
        self.l = n - Ls
        self.radix = len(alphabet)
        self.pblen = int(np.ceil(np.log(self.l) / np.log(self.radix))) #code length for positions
        self.lblen = int(np.ceil(np.log(self.Ls) / np.log(self.radix))) #code length for lengths
        
    """Pads the string from the left with a repetition
    of padsize times the symbol symb, where padsize is
    the length of the search buffer."""
    def iniPad(self, symb, padsize, cad):
        pad = padsize * symb
        return pad + cad

    """Looks for the longest string that starts at the 
    search buffer and that is a prefix of the lookahead
    buffer, using a SUFFIX TREE for the search. This method
    traverse the SuffixTree exploting the propertie of
    McCreights' SuffixTree construction: only the first symbol
    of the edge label must be checked in every node."""
    def ST_based_reproducible_extension(self, STD, s, lim=None, limsup=None):
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

    """Looks for the longest string that starts at the 
    search buffer and that is a prefix of the lookahead
    buffer."""
    def reproducible_extension(self, search, lookahead):
        pos = -1
        size = 0
        char = ''
        for prefixsize in range(1,min(self.l,len(lookahead))):
            prefix = lookahead[:prefixsize]
            p = search.find(prefix,0,self.l+prefixsize-1)
            #print(p, prefix, search[0:l+prefixsize-1])
            if p >= 0:
                pos = p
                size = prefixsize
                char = lookahead[size]
            else:
                break
        return pos, size, char


    """Obtain the code word associated with the prefix of
    length size, located at pos and followed by the character 
    char."""
    def codeWord(self, pos, size, char):
        #print(pos, size, char)
        cd = dec2radix(self.radix, pos).rjust(self.pblen, '0')
        cd += dec2radix(self.radix, size).rjust(self.lblen, '0')
        cd += char
        return cd

    """Code the string cad composed from symbols in alphabet
    using awindow of length n with a lookahead buffer of 
    length Ls."""
    def code(self, cad):    
        totcode = []
        pcad = self.iniPad(self.alph[0], self.l, cad)
        #print(pcad)
        wpos = 0
        while wpos < len(pcad)-self.Ls: #window position (does it need to go to the end?)
            pos, size, char = self.reproducible_extension(pcad[wpos:wpos+self.n], pcad[wpos+self.l:wpos+self.n])
            #print(pcad[wpos:wpos+self.n], " - ", pcad[wpos+self.l:wpos+self.n])
            #print(pos, size, char)
            c = ''
            if pos >= 0:
                c = self.codeWord(pos, size, char)
                #print('Found: ',c)
            else:
                c = self.codeWord(0, 0, pcad[wpos+self.l])
                #print("Prefix not found: ",c)
            
            totcode.append(c)
            wpos += max(1,size+1)
        return totcode

    """Code the string cad composed from symbols in alphabet
    using awindow of length n with a lookahead buffer of 
    length Ls, using the procedure specified in article:
    'Linear Algorithm for Data Compression via String Matching'
    MICHAEL RODEH, VAUGHAN R. PRATT, SHIMON EVEN. (Sec 4.1, pag 22)"""
    def ST_code(self, sigma, n, Ls):
        sigma += 'Z'
        
        totcode = []
        sigma = self.iniPad(self.alph[0], self.l-1, sigma)
        
        F = Ls
        N = n - Ls
        print("F:= ", F, " N:= ", N, ".")
        if not F < N/2:
            raise ValueError("Please make sure F < N/2")

        i = 0
        sigma_i = sigma[i*(N-F) : i*(N-F)+(N-1)]

        window_inf = 0
        window_sup = N-1
        while len(sigma[window_sup : window_sup + F]) > 0:
            #print("VENTANA:= ", sigma[window_inf : window_sup], " - ", sigma[window_sup : window_sup + F])

            if i*(N-F)+(N-1) < window_sup:
                i = i + 1
                sigma_i = sigma[i*(N-F) : i*(N-F)+(N-1)]
            elif (i+1)*(N-F) <= window_sup and window_sup <= i*(N-F)+(N-1):
                S_1i = STree(sigma[(i-1)*(N-F) : (i-1)*(N-F)+(N-1)])
                S_1i.root.traverse_postorder(postorder_scan_set)
                S_i = STree(sigma_i)
                S_i1 = STree(sigma[(i+1)*(N-F) : (i+1)*(N-F)+(N-1)])

                l_1i = self.ST_based_reproducible_extension(S_1i, sigma[window_sup : window_sup + F], lim=window_inf-(i-1)*(N-F))

                if l_1i is not None:
                    l_1i = (l_1i[0] - (window_inf-(i-1)*(N-F)), l_1i[1])

                l_i = self.ST_based_reproducible_extension(S_i, sigma[window_sup : window_sup + F], limsup=window_sup-i*(N-F))

                l_i = (l_i[0] + (i*(N-F)-window_inf), l_i[1])
                
                l_i1 = self.ST_based_reproducible_extension(S_i1, sigma[window_sup : window_sup + F], limsup=window_sup-(i+1)*(N-F))

                l_i1 = (l_i1[0] + ((i+1)*(N-F)-window_inf), l_i1[1])

                LIST = [(0,0), l_1i, l_i, l_i1]
                l = sorted([i for i in LIST if i is not None and i[0] < N-1], key=lambda x:x[1])[-1]

                #print("LIST:= ", LIST, ". from which we choosed: ", l)
                if l[1] == 0:
                    c = self.codeWord(0, 0, sigma[window_sup + l[1]])
                    totcode.append(c)
                    l = (0, 0)
                else:
                    #print(l, sigma[window_sup + l[1]])

                    c = self.codeWord(l[0], l[1], sigma[window_sup + l[1]])
                    totcode.append(c)

                window_inf = window_inf + l[1] + 1
                window_sup = window_sup + l[1] + 1
            else:
                S_1i = STree(sigma[(i-1)*(N-F) : (i-1)*(N-F)+(N-1)])
                S_1i.root.traverse_postorder(postorder_scan_set)
                S_i = STree(sigma_i)
                
                l_1i = self.ST_based_reproducible_extension(S_1i, sigma[window_sup : window_sup + F], lim=window_inf-(i-1)*(N-F))

                if l_1i is not None:
                    l_1i = (l_1i[0] - (window_inf-(i-1)*(N-F)), l_1i[1])

                l_i = self.ST_based_reproducible_extension(S_i, sigma[window_sup : window_sup + F], limsup=window_sup-i*(N-F))

                l_i = (l_i[0] + (i*(N-F)-window_inf), l_i[1])

                LIST = [(0,0), l_1i, l_i]
                l = sorted([i for i in LIST if i is not None and i[0] < N-1], key=lambda x:x[1])[-1]
                
                #print("LIST:= ", LIST, ". from which we choosed: ", l)
                if l[1] == 0:
                    c = self.codeWord(0, 0, sigma[window_sup + l[1]])
                    totcode.append(c)
                    l = (0, 0)
                else:
                    #print(l, sigma[window_sup + l[1]])

                    c = self.codeWord(l[0], l[1], sigma[window_sup + l[1]])
                    totcode.append(c)    

                window_inf = window_inf + l[1] + 1
                window_sup = window_sup + l[1] + 1

        return totcode

#Decoder-----------------------------------------------------------------
class Decoder:
    def __init__(self, alphabet, n, Ls):
        self.alph = alphabet
        self.n = n
        self.Ls = Ls
        self.l = n - Ls
        self.radix = len(alphabet)
        self.pblen = int(np.ceil(np.log(self.l) / np.log(self.radix))) #code length for positions
        self.lblen = int(np.ceil(np.log(self.Ls) / np.log(self.radix))) #code length for lengths
    
    """Obtain pos, size and char from a code."""
    def decodeWord(self, word):
        pos = radix2dec(self.radix, word[0:self.pblen])
        size = radix2dec(self.radix, word[self.pblen:self.pblen+self.lblen])
        char = word[self.pblen+self.lblen]
        return pos, size, char

    """Generate the current word from the decoded word."""
    def genWord(self, pcode, ibuff, cad):
        pos, size, char = self.decodeWord(pcode)
        #print(pos, size, char)
        ini = ibuff-self.l+pos
        sstr = ''
        i = 0
        while i < size:
            sstr += cad[ini + i%(ibuff-ini)] #If the word exceeds the buffer, complete by repeating characters from the begining of the prefix
            i += 1
        return sstr+char

    """Generate the current word from the decoded word."""
    def ST_genWord(self, pcode, ibuff, cad):
        pos, size, char = self.decodeWord(pcode)
        #print(pos, size, char)
        ini = ibuff-self.l+pos+1 # <------------------------------------------- +1 here
        sstr = ''
        i = 0
        while i < size:
            sstr += cad[ini + i%(ibuff-ini)] #If the word exceeds the buffer, complete by repeating characters from the begining of the prefix
            i += 1
        return sstr+char
    
    """Obtain the original string from the code."""
    def decode(self, totcode, ST=False):
        cad = self.l*self.alph[0]
        i = len(cad) #initial position of lookahead buffer
        for c in totcode:
            if ST:
                sstr = self.ST_genWord(c, i, cad)
            else:
                sstr = self.genWord(c, i, cad)
            #print(sstr)
            cad += sstr
            #print(cad)
            i += len(sstr)
        return cad[self.l:] #Remove the initial search buffer


#Functions---------------------------------------------------------------
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

"""Function for calculating the decimal representation
of rnum which is in base radix."""
def radix2dec(radix, rnum):
    dec = 0
    for i in range(0,len(rnum)):
        dec += int(rnum[len(rnum)-i-1])*(radix**i)
    return dec

def postorder_scan_set(self):
    if self.is_leaf():
        self.new_idx = self.idx
        return self.idx

    else:
        idxs = []
        for node in self.transition_links.values():
            idxs.append(postorder_scan_set(node))

        self.new_idx = max(idxs)
        return self.new_idx






alphabet = 'ATGC' # symbols in the alphabet

def codeWord(pos, size, n, Ls):
    radix = len(alphabet)
    pblen = int(np.ceil(np.log(n-Ls) / np.log(radix))) #code length for positions
    lblen = int(np.ceil(np.log(Ls) / np.log(radix))) #code length for lengths

    #print(pos, size, char)
    cd = dec2radix(radix, pos).rjust(pblen, '0')
    cd += dec2radix(radix, size).rjust(lblen, '0')
    #cd += char
    return cd

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

def codeString(string, n, Ls):
    totcode = []
    
    window_inf = 0
    window_sup = n
    window = string[window_inf:window_sup]
    while len(window) > 0:
        print(window)
        print(window[n-Ls:])
        ST = STree(window)
        pos, size = ST_based_reproducible_extension(ST, window[n-Ls:], limsup=n-Ls)

        print("pos=", pos, "long=", size)

        #c = codeWord(pos, long, n, Ls)
        totcode.append((pos, size, string[window_sup+size-1:window_sup+size]))

        if size == 0:
            size = 1

        window_inf += size
        window_sup += size

        window = string[window_inf:window_sup]
        
    print(totcode)


def list_to_string(L):
    S = ""
    for l in L:
        S += str(l)
    return S

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
    return SR

def encode_7_38(n):
    radix = 2
    decnum = n-7
    return dec2radix(radix, decnum)

def fibonacci(n):
    a, b = 0, 1
    for i in range(0, n):
        a, b = b, a + b
    return a

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
        cod_x = dec2radix(2, x)
    elif 7 <= x and x <= 38:
        cod_x = encode_7_38(x)
    else:
        cod_x = list_to_string(fibonacci_encode_number(x))

    if y < 7:
        cod_y = dec2radix(2, y)
    elif 7 <= y and y <= 38:
        cod_y = encode_7_38(y)
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
            print(S[i])
            encoded_symbol = encode_symbol(S[i])
            print("encoded symbol: ", encoded_symbol)
            totcode.append(encoded_symbol)
            size = 1
        else:
            print((pos, size))
            encoded_pair = encode_pair((pos, size))
            print("encoded pair: ", encoded_pair)
            totcode.append(encoded_pair)
        i = i+size

    print(totcode)
        
#---------------------------------------------------------------------------BEGGINING OF EXECUTION---------------------------------------------------------------------

cad = 'AACTGTTGTTGTTAGAACTGTTGTT'
print(len(cad))
"""
alphabet = 'ATGC' # symbols in the alphabet
n = 20 # window size
Ls = 5 # n/2 <= lookahead buffer size < n

coder = Coder(alphabet, n, Ls)
decoder = Decoder(alphabet, n, Ls)

c = coder.code(cad)
print(c)

ccad = decoder.decode(c)
print(ccad)

ST_c = coder.ST_code(cad, n, Ls)
print(ST_c)

ST_ccad = decoder.decode(ST_c, ST=True)
ST_ccad = ST_ccad[0:len(ST_ccad)-1]

print(cad)
print(ST_ccad)

print(ST_ccad == cad)
print(len(ST_c[0])*len(ST_c))
"""
#codeString(cad, 10, 5)
BioCompress(cad, cad)
