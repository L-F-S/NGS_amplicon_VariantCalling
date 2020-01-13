#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:40:49 2019

@author: L-F-S
@ University of Trento, Italy
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dataname="/home/lorenzo/cereseto lab/Aiello_EvolvR_Jan2020/03_variant_calling/L-JON_S2_L001_R1_001_trimmed.fq.gz.sorted.pileup"
pileup=pd.read_csv(dataname, sep="\t", header=0)
print(pileup.columns) #chr, pos, ref , A, C , G , T, af , cov
print(pileup.head(10))  #0-based
#%%
#index is the 1 based position (position 1 is missing for some reason, ma 
# fortunatamente non ci interessa)

#slice to only keep CJ_
CJ_start=151
CJ_end=3213

CJ_pileup=pileup[pileup["pos"]>=CJ_start]
CJ_pileup=CJ_pileup[CJ_pileup["pos"]<=CJ_end]
#rewrite positions:
CJ_pileup["pos"]=CJ_pileup["pos"]-CJ_start+1
print(CJ_pileup.tail)
#%%
#we are interested in values relative to coverage:
for letter in ["A", "T","C","G"]:
    pileup[letter]=round(pileup[letter]/pileup["cov"], 3)
print(pileup.tail())
#%%
print(pileup[pileup["pos"]==201])
#print(CJ_pileup["A"]+CJ_pileup["T"]+CJ_pileup["C"]+CJ_pileup["G"]) era un check: funge
#%%
#trucchetto per evidenziare il label:
# mi sevìrve per plottare più agilmente quello che voglio
def valore_ref(row): #1 aggiungi una colonna colla frazione della base reference
    letter=row["ref"]
    return row[letter]
def remove_ref(row):  #2 azzera la base reference
    letter=row["ref"]
    row[letter]=0.0
    return row

pileup["ref_value"]=pileup.apply(valore_ref, axis=1)
pileup=pileup.apply(remove_ref, axis=1)
#%%
pileup.tail()
#%%
#plot senza distanze vere
cutoff=0.0 #allelic fraction cutoff over which to consider a mutation
sliced=pileup[pileup["af"]>=cutoff]

#%%
to_print=sliced.drop(['chr', "ref", "af", "cov"], axis=1)

color_of={"T":"#0066ff", "A":"#ffd11a","C":"red", "G":"#00e600", "ref_value":"grey"}

plot_of={}
plt.figure(figsize=(101,20))
width=10
#per farlo con distanze giuste lungo il coso, sostituisci 'np.arange(len(to_print["pos"]))' con to_print["pos"], e mettiwidth=10 invece che 0.3
plot_of["ref_value"]=plt.bar(np.arange(len(to_print["pos"])), \
       to_print["ref_value"], width, color=color_of["ref_value"],)
previous_letter=to_print["ref_value"]
for letter in ["A","T","C","G"]:
    plot_of[letter]=plt.bar(np.arange(len(to_print["pos"])), to_print[letter], width, bottom=previous_letter, color=color_of[letter])
    previous_letter+=to_print[letter]
#plt.legend((plot_of["A"][0],plot_of["T"][0],plot_of["C"][0],plot_of["G"][0]),("A","T","C","G"), fontsize=50)
plt.yticks(size=50)

plt.xticks(np.arange(len(to_print["pos"])), [str(bu["pos"])+":"+str(bu["ref"]) for bo, bu in sliced.iterrows()], size=50)
plt.savefig("mutated_positions_all.png") 
#%%
#plot con distanze ingrandite maggiori sono le allelic fraction not yet
#%%
#translate to protein
#hm..interesting... I should call mutations first, i suppose...