# -*- coding: cp1250 -*-
from Bio import SeqIO
from Bio import Entrez
import Bio.Seq
from math import log
import random
import cPickle
import matplotlib.pyplot as plt
import pylab

#==========================================================================
def read(ime):
    """Iz datodeke ime preberi podake in jih vrni kot seznam."""
    f = open(ime)
    fc = f.read()
    return fc

#preberi podatke
data1=read("dnaC00C39.txt")
data2=read("dnaC40C79.txt")

#print seq

def read(ime, m):
    """Iz datodeke ime preberi podake in jih vrni kot seznam."""
    slovar={}
    matrika=[]
    f = open(ime)
    fc = f.read()
    fc = fc.split("\n")
    if m == 'd':
        for r in fc:
            s=r.split("\t")
            slovar[s[0]]=s[1]
        return slovar
    else:
        for r in fc:
            s=r.split(",")
            matrika.append(s)
        return matrika
        

data=read("aminoacids0-39.txt",'d')
amino=read("amino_file.txt",'a')


#==========================================================================
###ORF
##def codon_walk(s, frame=0):
##    for ix in range(frame, len(s), 3):
##        yield ix, s[ix:ix+3]
##
##
##def ORF(s,start_codon,stop_codon):             
##    genes=[]
##    for frame in range(3):
##        startP=None
##        codon_f=codon_walk(s, frame)
##        for st, codon in codon_f:
##            if codon in start_codon and startP==None:
##                startP=st
##            if codon in stop_codon and startP!=None:
##                genes.append((startP,st))
##                startP=None
##    return genes
##
##geni1=ORF(data1,amino[0][1:],amino[-1][1:])
##geni2=ORF(data1[::-1],amino[0][1:],amino[-1][1:])
##
##geni=geni1+geni2
####for g in geni:
####    print g
####    
##print len(geni)
##a= str(Bio.Seq.Seq(data1[23923:24202]).translate(table=6))
##b= data['C01']
##print a,'\n',b

#==========================================================================
##def prevedi(geni,data):
##    seznam=[]
##    for g in geni:
##        gen=data[g[0]:g[1]]
##        seznam.append(str(Bio.Seq.Seq(gen).translate(table=6)))
##    return seznam

def prevedi(data):
    seznam=[]
    for i in range(3):
        seznam.append(str(Bio.Seq.Seq(data[i:]).translate(table=6)))
    return seznam
        
seznam=prevedi(data1)+prevedi(data1[::-1])
print len(seznam)
##seznam = prevedi(geni1,data1)+prevedi(geni2,data1[::-1])
##print seznam[0]
#==========================================================================

def read_table(fname):
    d = {}
    # key: (a1, a2), value: "penalty" score
    lines = open(fname, "rt").readlines()
    alpha = lines[0].rstrip('\n\r').split()
    assert(len(alpha) == len(lines)-1)
    for r in lines[1:]:
        r = r.rstrip('\n\r').split()
        a1 = r[0]
        for a2, score in zip(alpha, r[1:]):
            d[(a1, a2)] = int(score)
    return d

blosum50 = read_table("blosum50.txt")
delta_book = {}
for a1, a2 in blosum50.keys():
    if a1 == a2:
        delta_book[(a1, a2)] = 1
    else:
        delta_book[(a1, a2)] = -1

def align_nw(s, t, delta, gap_p=-1):
    M = [[0]*(len(t)+1) for i in range(len(s)+1)]
    pred = [[(0,0)]*(len(t)+1) for i in range(len(s)+1)]
    for i in range(1, len(s)+1):
        M[i][0] = M[i-1][0] + gap_p
        pred[i][0] = (i-1, 0)

    for j in range(1, len(t)+1):
        M[0][j] = M[0][j-1] + gap_p
        pred[0][j] = (0, j-1)

    for i in range(1, len(s)+1):
        for j in range(1, len(t)+1):
            M[i][j], pred[i][j] = \
                max((M[i-1][j] + gap_p, (i-1, j)),
                (M[i][j-1] + gap_p, (i, j-1)),
                (M[i-1][j-1] + delta[(s[i-1], t[j-1])], (i-1, j-1)))
    return M, pred

#=============================================================================
def align_sw(s, t, delta, gap_p=-1):
    Ma = [[0]*(len(t)+1) for i in range(len(s)+1)]
    pred = [[(0,0)]*(len(t)+1) for i in range(len(s)+1)]

    for i in range(1, len(s)+1):
        for j in range(1, len(t)+1):
            Ma[i][j], pred[i][j] = \
                max((Ma[i-1][j] + gap_p, (i-1, j)),
                (Ma[i][j-1] + gap_p, (i, j-1)),
                (Ma[i-1][j-1] + delta[(s[i-1], t[j-1])], (i-1, j-1)),
                (0, (0,0)))
    return Ma, pred

def traceback_nw(s, t, pred):
    walk = [(len(s), len(t))]
    prev_i, prev_j = walk[-1]
    while pred[prev_i][prev_j] != (0, 0):
        walk.append(pred[prev_i][prev_j])
        prev_i, prev_j = walk[-1]
    walk.reverse()
    return walk

def traceback_sw(s, t, pred, mat, zac_celica=None):
    walk = [zac_celica]
    prev_i, prev_j = walk[-1]
    while pred[prev_i][prev_j] != (0, 0):
        prev_i, prev_j = pred[prev_i][prev_j]
        if mat[prev_i][prev_j] == 0:
            break
        walk.append((prev_i, prev_j))
    walk.reverse()
    return walk

def pp_alignment(s, t, walk):
    def pp(indx, n):
        ret = []
        last = None
        for i in indx:
            if not i or (i == last):
                ret.append('-')
            else:
                ret.append(n[i-1])
            last = i
        return ret
    ss = [i for i, _ in walk]
    tt = [j for _, j in walk]
    ps = pp(ss, s)
    pt = pp(tt, t)

    ps = "".join(ps)
    pt = "".join(pt)
    print ps
    print "".join([' ', '|'][ns == nt] for ns, nt in zip(ps, pt))
    print pt
    
#=============================================================================
def racunaj(s,t):
##    s = 'VIVALASVEGAS'
##    t = 'VIVADAVIS'
##    s=seq
##    t=data["C01"]
##    print 's:', s,
##    print 't:', t
##    print
    import time
##    t1 = time.time()
##    mat, pr = align_nw(s, t, blosum50)
##    print "TIME", time.time() - t1
##
##    w = traceback_nw(s, t, pr)
##    print "score of the global alignment:", mat[-1][-1]
##
##    pp_alignment(s, t, w)
##    print

#________________________________________________
    print "local alignment"
    mat, pr = align_sw(s, t, delta_book) #blosum50)
    print len(mat), mat[30]
    loc_score = max(max(r) for r in mat)
    print "score (of the best) local alignment:", loc_score

    print

    
    for i, r in enumerate(mat):
        if loc_score in r:
            j = r.index(loc_score)
            print "possible local alignment with score:", mat[i][j]
            w = traceback_sw(s, t, pr, mat, (i, j))
            pp_alignment(s, t, w)
            print

        

##for i in range(len(seq)):
##        racunaj(seq[i])
dobri=["C01"]#,"C03","C05","C08","C25","C36","C29"]
for d in dobri:#data.keys():
    for s in seznam:
        print d, len(s), len(data[d]) 
        racunaj(s, data[d])

##        for i in range(len(seq_r)):
##            racunaj(seq_r[i],data[d])


##racunaj(b,a)    
