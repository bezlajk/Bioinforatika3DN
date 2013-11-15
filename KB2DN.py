# -*- coding: cp1250 -*-
from Bio import SeqIO
from math import log
from Bio import Entrez
import random
import cPickle
import matplotlib.pyplot as plt
import pylab

#id = "NC_007346"
id = "NC_006058"
loaded = True
#==========================================================================
#Pridobi sekvenco
if not loaded:
    handle = Entrez.efetch(db="nucleotide", rettype="gb", id=id,
                           email="a@uni-lj.si")
    rec = SeqIO.read(handle, "gb")
    handle.close()
    cPickle.dump(rec, file(id, "w"))
else:
    rec = cPickle.load(file(id))
    
data = str(rec.seq)
data_rc = str(rec.seq.reverse_complement())
#==========================================================================
print "Sequence length", len(data)

def codon_walk(s, frame=0):
    for ix in range(frame, len(s), 3):
        yield ix, s[ix:ix+3]

#===========================================================================
#verjetnost stopa
def probabilityStop(data,data_rc):
    n_stop = 0
    for s in (data, data_rc):
        for frame in range(3):
            n_stop += sum(1 for _, codon in codon_walk(s, frame) if codon in stop_codon)
    return float(n_stop) / (2 * (len(data) - 2))
#==========================================================================
#ORF
def ORF(s,start_codon,stop_codon):             
    genes=[]
    for frame in range(3):
        startP=None
        codon_f=codon_walk(s, frame)
        for st, codon in codon_f:
            if codon in start_codon and startP==None:
                startP=st
            if codon in stop_codon and startP!=None:
                genes.append((startP,st))
                startP=None
    return genes

#===============================================================================
# Real ORF
lengths = [(f.location.end.position-f.location.start.position)/3-1 for f in rec.features
           if f.type=="CDS"]

#===============================================================================
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
        if seq[i:i+3]==start_codon and len(gen_start)<4 and intron==[]:
            #našli smo start kodon, lahko je prvi ali pa tudi ne
            if gen_start==[]:
                gen_start.append((i%3,i))
            else:
                frame=[f for f,_ in gen_start]
                if i%3 not in frame:
                    gen_start.append((i%3,i))

        if seq[i:i+2]==start_i and len(istart)<4:
            #našli smo kandidata za intron
            if gen_start!=[] and istart==[]:
                istart.append((i%3,i))
            
            elif gen_start!=[] and istart!=[]:  
                frame_i=[f for f,_ in gen_start]
                if i%3 not in frame_i:
                    istart.append((i%3,i))
                       
        if seq[i:i+2]==stop_i:
            #našli smo konec introna
            [intron.append((j,i)) for _,j in istart ]
            istart=[]
                
        if seq[i:i+3]==stop_codon and gen_start!=[] and istart==[]:
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

#===============================================================================
def Matching(trashold,orf,orf_rc):
    CDS=[(f.location.start.position,f.location.end.position) for f in rec.features if f.type=="CDS" and f.strand==1]
    CDS_rc=[(f.location.end.position,f.location.start.position) for f in rec.features if f.type=="CDS" and f.strand==-1]
    CDS_start=[f.location.start.position for f in rec.features if f.type=="CDS" and f.strand==1]
    CDS_start_rc=[f.location.end.position for f in rec.features if f.type=="CDS" and f.strand==-1]
    orf_start=[b for b,e in orf if (e-b)/3>=trashold]
    orf_start_rc=[len(data)-b for b,e in orf_rc if (e-b)/3>=trashold]
    
    matches_s=len(list(set(CDS_start)&set(orf_start)))
    matches_rc_s=len(list(set(CDS_start_rc)&set(orf_start_rc)))
    
    matches=len(list(set(CDS)&set(orf_start)))
    matches_rc=len(list(set(CDS_rc)&set(orf_start_rc)))
    
    dr=len(CDS)+len(CDS_rc)
    dp=len(orf_start)+len(orf_start_rc)
    
    precision1=(matches_s+matches_rc_s)/float(dp)
    recall1=(matches_s+matches_rc_s)/float(dr)

    precision2=(matches+matches_rc)/float(dp)
    recall2=(matches+matches_rc)/float(dr)

    return precision1,recall1,precision2,recall2
#=======================================================================================================
#=======================================================================================================

if id =="NC_006058":
    start_codon="ATG"
    stop_codon = "TGA" # in M genitalium TGA is not a stop codon
    
    #==================================
    p_stop=probabilityStop(data,data_rc)
    alpha = 0.010
    k = log(alpha) / log(1 - p_stop) + 1 # +1 for inclusion of stop codon
    print "Stop codon probability: %6.4f" % p_stop
    print "Minimal ORF length (alpha=%5.3f): %d" % (alpha, k)
    
    #==================================
    orf=ORF(data,start_codon,stop_codon)
    orf_rc=ORF(data_rc,start_codon,stop_codon)
    print "========================================================="
    print "ORF positive: ",len(orf)
    print "ORF negative: ",len(orf_rc)
    print "All ORF: ",len(orf)+len(orf_rc)
    #==================================
    
    lengths = [(stop-start)/3 for (start,stop) in orf+orf_rc]
    plt.figure("ORFS Paramecium tetraurelia")
    plt.title("ORFS Paramecium tetraurelia")
    plt.hist([l for l in lengths if l>100], 40)
    
    #==================================
    print "========================================================="
    print "Total features: %d"% len([a for a in rec.features if f.type=="CDS"])
    plt.figure("CDR Paramecium tetraurelia")
    plt.title("CDS Paramecium tetraurelia")
    plt.hist(lengths, 40)
    
    #==================================
    neki=findIntron(data,start_i,stop_i,start_codon,stop_codon)
    neki1=findIntron(data_rc,start_i,stop_i,start_codon,stop_codon)
    print "============================================================"
    geni=set(neki)|set(neki1)
    print len(geni),"geni"

    #=================================
    trashold=125
    P1=Matching(trashold,orf,orf_rc)

    prec=[]
    reca=[]
    gap=10
    stop=800

    for i in range(0,stop,gap):
        p = Matching(i, orf, orf_rc)
        prec.append(p[0])
        reca.append(p[1])


    plt.figure("Precison and recall")
    plt.plot(range(0,stop,gap),prec)
    plt.plot(range(0,stop,gap),reca)
    plt.legend(["precision","recall"])

    print "========================================================="
    print "Precision and recall for ORF for trashold %d: "%trashold+str(P1[0:2])
    print "Precision and recall for ORF for trashold %d start and stop: "%trashold+str(P1[2:])
    
    #=================================
    prec=[]
    reca=[]
    gap=10
    stop=800

    for i in range(0,stop,gap):
        p = Matching(i, neki, neki1)
        prec.append(p[0])
        reca.append(p[1])


    plt.figure("Precison and recall modify")
    plt.plot(range(0,stop,gap),prec)
    plt.plot(range(0,stop,gap),reca)
    plt.legend(["precision","recall"])

    P2=Matching(trashold,neki,neki1)
    print "Precision and recall for ORF for trashold %d: "%trashold+str(P2[0:2])
    print "Precision and recall for ORF for trashold %d start and stop: "%trashold+str(P1[2:])
    pylab.show()

else:
    start_codon=["ATG", "TTG", "CTG"]
    stop_codon=["TAA", "TAG", "TGA"]
    
    #==================================
    p_stop=probabilityStop(data,data_rc)
    alpha = 0.010
    k = log(alpha) / log(1 - p_stop) + 1 # +1 for inclusion of stop codon
    print "Stop codon probability: %6.4f" % p_stop
    print "Minimal ORF length (alpha=%5.3f): %d" % (alpha, k)
    
    #==================================
    orf=ORF(data,start_codon,stop_codon)
    orf_rc=ORF(data_rc,start_codon,stop_codon)
    print "========================================================="
    print "ORF positive: ",len(orf)
    print "ORF negative: ",len(orf_rc)
    print "All ORF: ",len(orf)+len(orf_rc)

    #==================================
    lengths = [(stop-start)/3 for (start,stop) in orf+orf_rc]
    plt.figure("ORFS Emiliania huxleyi virus 86")
    plt.title("ORFS Emiliania huxleyi virus 86")
    plt.hist([l for l in lengths if l>100], 40)
    
    #==================================
    #print "========================================================="
    #print "Total nunber od gens: %d"% len([a for a in rec.features if f.type=="CDS"])
    plt.figure("CDR Emiliania huxleyi virus 86")
    plt.title("CDS Emiliania huxleyi virus 86")
    plt.hist(lengths, 40)
    #==================================

    trashold=125
    P1=Matching(trashold,orf,orf_rc)
    
    prec=[]
    reca=[]
    gap=10
    stop=400

    for i in range(0,stop,gap):
        p = Matching(i, orf, orf_rc)
        prec.append(p[0])
        reca.append(p[1])

    plt.figure("Precison and recall")
    plt.plot(range(0,stop,gap),prec)
    plt.plot(range(0,stop,gap),reca)
    plt.legend(["precision","recall"])

    print "========================================================="
    print "Precision and recall for ORF for trashold %d: "%trashold+str(P1[0:2])
    print "Precision and recall for ORF for trashold %d start and stop: "%trashold+str(P1[2:])

pylab.show()    
            
 
