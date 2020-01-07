#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 22:00:55 2017
@author: zw5j
"""

import os,sys,re,argparse,bisect
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=14
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.6)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")

def plot_occurrene_vs_motif(df):
    #.figure(figsize=(5,5))
    fig,ax1 = plt.subplots(figsize=(3.6,3.2))
    a1 = ax1.plot(df.index,df['#bindings'],label='bindings',color='navy')
    ax1.set_yscale('log')
    ax1.set_ylabel('# CTCF binding sites',fontsize=18)
    ax1.set_xlabel('Occupancy score',fontsize=18)
    ax1.set_xticks([0,250,500,750])
    #ax1.legend(bbox_to_anchor=(1.2,1))

    ax2 = ax1.twinx()
    a2 = ax2.plot(df.index,100*df['%motif'],label='% motif',color = 'coral')
    #ax2.legend(bbox_to_anchor=(1.2,0.7))

    ls = a1+a2
    labs = [l.get_label() for l in ls]
    ax1.legend(ls,labs,frameon=False,bbox_to_anchor=(.7,.25),fontsize=16,borderaxespad=.2,labelspacing=.2,handletextpad=0.5,handlelength=1.5)
    ax1.tick_params(axis='y',direction='out', length=3, width=.8, colors='black')
    ax1.tick_params(axis='x',direction='out', length=3, width=.8, colors='black')
    ax2.tick_params(axis='y',direction='out', length=3, width=.8, colors='black')
    ax2.set_ylabel('% bindings with motif',fontsize=18)
    plt.savefig('figures/unionCTCF_occupancy_score_vs_motif.pdf',bbox_inches = 'tight',transparent=True)
    plt.close()

    

def main():
    
    inf='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/binding_occupancy_with_motif.csv'
    out_df = pd.read_csv(inf,index_col=0)
    plot_occurrene_vs_motif(out_df)
    



if __name__ == '__main__':
  
    main()

