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
#data1=read("dnaC00C39.txt")
data1=read("dnaC40C79.txt")

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


data=read("aminoacids40-79.txt",'d')
amino=read("amino_file.txt",'a')
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

data1_rc=str(Bio.Seq.Seq(data1).reverse_complement())
seznam=prevedi(data1)+prevedi(data1_rc)
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
    return walk, prev_i+1

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

#==============================================================================
def racunaj_globalno(s,t,risi):
    import time
    t1 = time.time()
    mat, pr = align_nw(s, t, blosum50)
    #print "TIME", time.time() - t1

    w = traceback_nw(s, t, pr)
    if risi==1:
        print "score of the global alignment:", mat[-1][-1]
        pp_alignment(s, t, w)
    return mat[-1][-1]
#==============================================================================
def racunaj_lokalno(s,t,koef):
    #print "local alignment"
    if koef=='a':
        mat, pr = align_sw(s, t, delta_book)
    else:
        mat, pr = align_sw(s, t, blosum50,koef)
    loc_score = max(max(r) for r in mat)
    #print "score (of the best) local alignment:", loc_score
    return mat, pr, loc_score

def izpisi(s,t,mat,pr,loc_score,risi):
    min_dolzina=10000000
    for i, r in enumerate(mat):
        if loc_score in r:
            j = r.index(loc_score)
            w, z = traceback_sw(s, t, pr, mat, (i, j))
            #print "possible local alignment with score:", mat[i][j]
            #print "začetek sekvence je na mestu: ", z, "dolžina sekvence pa je: ", len(w), loc_score
            #if len(w)<min_dolzina:
            #    min_dolzina=len(w[0])
            #    pravi=[z,len(w[0]),w]
            #print z,w,w[0][1],len(w)
            wt=w[0][1]
    if risi==1:
        print "začetek sekvence je na mestu: ", z, "dolžina sekvence pa je: ", len(w)
        pp_alignment(s, t, w)#pravi[2])
    return z,len(w), wt #pravi[0], pravi[1]

#==============================================================================
def rekurzija(z, k, th, n, s, t):#th=trachhold, n natančnost
    #print z,konec
    if n!='a':
        n-=5
    if abs(z-k) < th: return []
    else:
        mat,pr,m = racunaj_lokalno(s[z:k], t, n)
        zacetek, dolzina,_ = izpisi(s[z:k], t, mat, pr, m, 0)
        #print zacetek, z, z+zacetek, dolzina, z+zacetek+dolzina, k
        ekson1 = rekurzija(z, z+zacetek, th, n, s, t)
        ekson2 = rekurzija(z+zacetek+dolzina, k, th, n, s, t)
        return ekson1+ekson2+[[z+zacetek, z+zacetek+dolzina]]
#==============================================================================

seznam_la=[]
dobri=data.keys()
#dobri=['C01']#,'C01','C02']
for d in dobri:
    #print '\n Analiza gena %s \n'%d
    maxi=None
    maxi_mat=[]
    maxi_pr=[]
    maxi_s=[]
    maxi_l=0
    for i, s in enumerate(seznam):
        #print d, i, len(data[d])
        mat,pr,m = racunaj_lokalno(s, data[d],-10)
        z, l, _= izpisi(s,data[d],mat,pr,m,0)
        if m>maxi or (maxi==m and l>maxi_l):
            maxi=m
            maxi_mat=mat
            maxi_pr=pr
            maxi_s=s
            maxi_i=i
            maxi_l=l
    z, l, zt = izpisi(maxi_s, data[d], maxi_mat, maxi_pr, maxi, 0)
    dolzina_gena=len(data[d])

    koef=-30
    zacetek=z-zt-5
    konec=z+dolzina_gena-zt
    #print s[8046:8219]
    #print z, z+zt, s[z:z+l]
    #ekson = rekurzija(zacetek, konec, dolzina_gena*0.20, -100, s, data[d])
    #print ekson
    #ekson=sorted(ekson)
    #ekstroni=[z-zt]
    #i=0
    #j=0
    #while i<(len(ekson)-1):
    #    if ekson[i][0]>ekstroni[0]:
    #        if abs(ekson[i][0]-ekson[i][1])>=3:
    #            if abs(ekson[j][1]-ekson[i+1][0])>10:
    #                ekstroni.append(ekson[i][1])
    #                ekstroni.append(ekson[i+1][0])
    #           j=i
     #       i+=1
    #ekstroni.append(ekson[-1][1])
    seznam_la.append([d,[z+zt,z+l+zt]])

for sez in seznam_la:
    print "%s"%sez[0],
    for i in range(0,len(sez[1])):
        print "\t%d"%(sez[1][i]*3+maxi_i),
    print "\n",
print