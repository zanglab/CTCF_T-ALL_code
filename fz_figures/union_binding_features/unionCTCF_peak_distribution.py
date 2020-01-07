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
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
import CTCF_TALL_modules_new 
matplotlib.rcParams["font.sans-serif"] = ["Arial"]




def plt_pdf_sns(mylist,bins=0,ylog=False,figname=None):
    
    plt.figure(figsize=(3,3))   
    #sns.distplot(mylist,kde=True,hist_kws = dict(color='b'))
    plt.hist(mylist,bins=bins,color='navy',rwidth=1.05)   
    plt.axes().tick_params(axis='y',direction='out', length=3, width=.8, colors='black')
    plt.axes().tick_params(axis='x',direction='out', length=3, width=.8, colors='black')
    #plt.title('Distribution of peak numbers',fontsize=20)
    plt.ylabel('Number of datasets',fontsize=16)
    plt.xlabel('Number of peaks',fontsize=16)
    plt.xlim([2000,150000])
    plt.xticks([2000,50000,100000,150000],rotation=0,ha='center') 
    if figname:
        plt.savefig(figname,bbox_inches='tight',pad_inches=0.02,transparent=True)   
    plt.close()   


    

collection_df = CTCF_TALL_modules_new.return_collection_df()
# print(collection_df)
figname ='figures/unionCTCF_peak_distribution.pdf'
plt_pdf_sns(collection_df['peak_nums'].values,bins=40,figname=figname)
