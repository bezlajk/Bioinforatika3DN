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
#ORF
def codon_walk(s, frame=0):
    for ix in range(frame, len(s), 3):
        yield ix, s[ix:ix+3]


def ORF(s,start_codon,stop_codon):             
    genes=[]
    geni=[]
    for frame in range(3):
        startP=None
        codon_f=codon_walk(s, frame)
        for st, codon in codon_f:
            if codon in start_codon and startP==None:
                startP=st
            if codon in stop_codon and startP!=None:
                genes.append((startP,st))
                startP=None
        geni.append(genes)
    return geni

geni1=ORF(data1,amino[0][1:],amino[-1][1:])
geni2=ORF(data1[::-1],amino[0][1:],amino[-1][1:])

geni=geni1+geni2

Geni=[]
for i in range(len(geni)):
    ge=[]
    for g in geni[i]:
        ge.append([g[0]/3,str(Bio.Seq.Seq(data1[g[0]:g[1]]).translate(table=6))])
    Geni.append(ge)
print Geni
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

def traceback_nw(s, t, pred):
    walk = [(len(s), len(t))]
    prev_i, prev_j = walk[-1]
    while pred[prev_i][prev_j] != (0, 0):
        walk.append(pred[prev_i][prev_j])
        prev_i, prev_j = walk[-1]
    walk.reverse()
    return walk

#=============================================================================
def align_sw(s, t, delta, gap_p=-1):
    M = [[0]*(len(t)+1) for i in range(len(s)+1)]
    pred = [[(0,0)]*(len(t)+1) for i in range(len(s)+1)]

    for i in range(1, len(s)+1):
        for j in range(1, len(t)+1):
            M[i][j], pred[i][j] = \
                max((M[i-1][j] + gap_p, (i-1, j)),
                (M[i][j-1] + gap_p, (i, j-1)),
                (M[i-1][j-1] + delta[(s[i-1], t[j-1])], (i-1, j-1)),
                (0, (0,0)))
    return M, pred


def traceback_sw(s, t, pred, mat, zac_celica=None):
    #tuki najdemo začetek niza
    walk = [zac_celica]
    prev_i, prev_j = walk[-1]
    while pred[prev_i][prev_j] != (0, 0):
        prev_i, prev_j = pred[prev_i][prev_j]
        if mat[prev_i][prev_j] == 0:
            break
        walk.append((prev_i, prev_j))
    walk.reverse()
    return walk, prev_i
#======================================================================
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
def racunaj_globalno(s,t,risi):
##    s = 'VIVALASVEGAS'
##    t = 'VIVADAVIS'
##    s=seq
##    t=data["C01"]
##    print 's:', s,
##    print 't:', t
##    print
    import time
    t1 = time.time()
    mat, pr = align_nw(s, t, blosum50)
    print "TIME", time.time() - t1

    w = traceback_nw(s, t, pr)
    print "score of the global alignment:", mat[-1][-1]
    if risi==1:
        pp_alignment(s, t, w)
    return mat[-1][-1]

#________________________________________________

def racunaj_lokalno(s,t):
    #print "local alignment"
    mat, pr = align_sw(s, t, delta_book) #blosum50)
    loc_score = max(max(r) for r in mat)
    #print "score (of the best) local alignment:", loc_score
    return mat, pr, loc_score

def izpisi(s,t,mat,pr,loc_score):
    for i, r in enumerate(mat):
        if loc_score in r:
            j = r.index(loc_score)
            w, z = traceback_sw(s, t, pr, mat, (i, j))
            print "possible local alignment with score:", mat[i][j]
    print "začetek sekvence je na mestu: ", z, "dolžina sekvence pa je: ", len(w)

    pp_alignment(s, t, w)
    return z, len(w)

#=========================================================================      

##def poisci_start(s,t,z,w,geni):
##    a=None
##    stev=0
##    while True:
##        if geni[stev][0]>z:
##            break
##        stev+=1
##    for g in geni:
##        print g
##    for g in geni:
##        mat2,pr2,m2=racunaj_lokalno(g[1],t)
##        izpisi(g[1],t,mat2,pr2,m2)
##    
def poisci_start(s,t,z,w):
    maxcena=None
    mesto=0
    mesto_o=0
    i=z
    j=z
    k=z+w
    while s[k]!='*':
        k+=1
    while True:
        j-=1
        if s[j]=='M':
            cena=racunaj_globalno(s[j:k],t,0)
            if maxcena<cena:
                print 'tuki\n'
                maxcena=cena
                mesto=j
            else: break
    print 'mesto začetka:', mesto, 'mesto konca:', k
    racunaj_globalno(s[mesto:k],t,1)
    while True:
        if i<k:
            i+=1
            if s[i]=='M':
                cena=racunaj_globalno(s[i:k],t,0)
                if maxcena<cena:
                    print 'tuki\n'
                    maxcena=cena
                    mesto_o=i
                else: break
        else:
            mesto_o=None
            break
    if mesto_o!=None:
        print 'mesto začetka:', mesto_o, 'mesto konca:', k
        racunaj_globalno(s[mesto_o:k],t,1)
             
    
    
    
#=========================================================================

print seznam[2][7974:8220]

dobri=["C01","C03","C05","C08","C25","C36","C29"]
for d in dobri:#data.keys():
    maxi=None
    maxi_mat=[]
    maxi_pr=[]
    maxi_s=[]
    for i, s in enumerate(seznam):
        #print d, i, len(data[d]) 
        mat,pr,m = racunaj_lokalno(s, data[d])
        if m>maxi:
            maxi=m
            maxi_mat=mat
            maxi_pr=pr
            maxi_s=s
            maxi_i=i
        
    z, w = izpisi(maxi_s,data[d],maxi_mat,maxi_pr,maxi)
    print
    poisci_start(maxi_s,data[d],z,w)#,Geni[maxi_i])



