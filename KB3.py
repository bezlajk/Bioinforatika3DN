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
start_i="GT"
stop_i="AG"

def kombiniraj(sez):
    sez1=[]
    for i,j in sez:
        for k,l in sez:
            if i<l:
                sez1.append((i,l))
            else:
                sez1.append((l,i))
    return sez1
            
def findIntron(seq,start_i,stop_i,start_codon,stop_codon):
    output=[]
    gen_start=[]
    intron=[]
    istart=[]
    stop=None
    for i,a in enumerate(seq):
        if seq[i:i+3] in start_codon and len(gen_start)<4 and intron==[]:
            #našli smo start kodon, lahko je prvi ali pa tudi ne
            if gen_start==[]:
                gen_start.append((i%3,i))
            else:
                frame=[f for f,_ in gen_start]
                if i%3 not in frame:
                    gen_start.append((i%3,i))

        if seq[i:i+2] in start_i and len(istart)<4:
            #našli smo kandidata za intron
            if gen_start!=[] and istart==[]:
                istart.append((i%3,i))
            
            elif gen_start!=[] and istart!=[]:  
                frame_i=[f for f,_ in gen_start]
                if i%3 not in frame_i:
                    istart.append((i%3,i))
                       
        if seq[i:i+2] in stop_i:
            #našli smo konec introna
            [intron.append((j,i)) for _,j in istart ]
            istart=[]
                
        if seq[i:i+3] in stop_codon and gen_start!=[] and istart==[]:
            if intron!=[]:
                introni=kombiniraj(intron)
                for _,l in gen_start:
                    for j,k in introni:
                        d=(j-l)+(i-k)
                        if d%3==0 and d<125*3:
                            output.append((l,i))
                            stop=1
                            break
                    if (l-i)%3==0 and stop!=1:
                        output.append((l,i))
                        stop=None
                    elif stop==1:
                        stop=None
                        break
                intron=[]
            else:
                for _,j in sorted(gen_start):
                    if (j-i)%3==0:
                        output.append((j,i))
                        gen_start=[]
                intron=[]
    return output

geni1=findIntron(data1,start_i,stop_i,amino[0][1:],amino[-1][1:])
geni2=findIntron(data1,start_i,stop_i,amino[0][1:],amino[-1][1:])
geni=geni1+geni2
print len(geni)
#==========================================================================
def pretvori(data):
    seq =[]
    for i in range(3):
        data1=data[i:]
        seq.append(str(Bio.Seq.Seq(data1).translate(table=6)))
    return seq

seq=pretvori(data1)
seq_r=pretvori(data1[::-1])



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

##    w = traceback_nw(s, t, pr)
##    print "score of the global alignment:", mat[-1][-1]

    #pp_alignment(s, t, w)
    print

#________________________________________________
    print "local alignment"
    mat, pr = align_sw(s, t, delta_book) #blosum50)
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
dobri=["C01","C03","C05","C08","C25","C36","C29"]
for d in dobri:#data.keys():
    for i in range(len(seq)):
        print d, i, len(data[d]) 
        racunaj(seq[i], data[d])

    for i in range(len(seq_r)):
        racunaj(seq_r[i],data[d])


    
