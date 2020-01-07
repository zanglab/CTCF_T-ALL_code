#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 20:12:42 2018

@author: zw5j
"""

import os,sys,argparse,glob,re,bisect
import numpy as np
import pandas as pd
from collections import Counter
import operator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#matplotlib.rcParams['agg.path.chunksize'] = 10000
matplotlib.rcParams['font.size']=20
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.6)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
import scipy.optimize
sns.set_style("ticks")



def fit_exp_nonlinear(t, y):
    opt_parms, parm_cov = scipy.optimize.curve_fit(model_function, t, y, maxfev=200000)
#    print(opt_parms, parm_cov)
    A,m, K= opt_parms
    return A,m, K

# def fit_exp_linear(t,y,C=0):
#     y = y-C
#     y = np.log(y)
#     K,log_A = np.polyfit(t,y,1)
#     A = np.exp(log_A)
#     return A,K,C

def model_function(t,A,m,K):
    #return A*np.exp(K*t)
#    return A*((t-m)**(K))
    return (1/A)*((t-m)**(K))

occupancy_count='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f2_occupancy_on_union_CTCFs/f2_occupancy_score/occupancy_count.csv'
df = pd.read_csv(occupancy_count,sep='\t',index_col=0)
x,y = df.index,df['#bindings'].values

a,b = 3,800
A,m,K = fit_exp_nonlinear(np.array(x[a:b]),np.array(y[a:b]));print('A,m,k\t',A,m,K)
fit_y = model_function(np.array(x),A,m,K)
#print(len(fit_y))
plt.figure(figsize=(3.8,3.5))
plt.plot(x,y,label='Observed',color='navy')
plt.plot(x,fit_y,label='Expected',color='limegreen')   
plt.axvline(x=771*.8,ymax=.3,ymin=0,color='grey',linestyle='--')
plt.legend(fontsize=16,borderaxespad=0.1,labelspacing=.1,handletextpad=0.2,loc="upper right",frameon=False)
plt.yscale('log')
# plt.xlabel('# Occupancy score\nin 771 datasets',fontsize=20)
plt.xlabel('Occupancy score',fontsize=20)
plt.ylabel('# CTCF binding sites',fontsize=20)
plt.axes().tick_params(axis='x',direction='out', length=4, width=1, colors='black')    
plt.axes().tick_params(axis='y',direction='out', length=4, width=1, colors='black')    
plt.savefig('figures/unionCTCF_occupancy_score_power_model_fit.pdf',bbox_inches = 'tight',pad_inches = 0.1,dpi=600,transparent=True)
plt.show()
plt.close()
sum_ori,sum_fit,flag = 0,0,0
for i in np.arange(len(x)-1,540,-1):
    yfit = model_function(i,A,m,K)
    sum_ori += y[i]
    sum_fit += yfit
#     flag = max(flag,sum_ori/sum_fit)
#    if flag>10: # sum_ori/sum_fit increase first then decrease
#        if sum_ori/sum_fit<8:
#     if i <480:
    print(i,y[i],yfit,(sum_ori)/sum_fit)
            #break

