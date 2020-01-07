import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
import re,bisect
import CTCF_TALL_modules
import scipy
from scipy import stats
sys.path.insert(0,os.path.abspath('modules'))
# import return_differential_binding
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
matplotlib.rcParams["font.sans-serif"] = ["Arial"]


def mark_pvalue(compr_pos,positions,box_vals,flag):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    print('\n',flag,'\n',s,p)
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),99)*1.05 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),0)*0.95
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    #if compr_pos[0]%3==0:
    #    y = y*1.5
    if p<0.05:
        if compr_pos[2] == 't':
            plt.plot([x1, x1, x2, x2], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*h, "{:.1e}".format(p), ha='center', va='bottom', color=col,fontsize=14)
        else:
            plt.plot([x1, x1, x2, x2], [y2, y2*1.1, y2*1.1, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*1.3, "{:.1e}".format(p), ha='center', va='bottom', color=col,fontsize=14)


def enrichment_compr_bar(baselists,targets,ylabel,ylim,xticklabels,colors,figname):
    
    plt.figure(figsize=(1.5,2.6))
    width = 0.6
    position=0
    a1 = len(baselists[0])
    a2 = len(set(baselists[0]).intersection(targets));print(ylabel);print(xticklabels[position],'\t',a1,'\t',a2)
    plt.bar(position,100*a2/a1,width=width,color = colors[position]);position+=1

    for baselist in baselists[1:]:
        b1 = len(baselist)
        b2 = len(set(baselist).intersection(targets));print(xticklabels[position],'\t',b1,'\t',b2)
        plt.bar(position,100*b2/b1,width=width,color = colors[position])
        s,p = stats.fisher_exact([[a1-a2,a2],[b1-b2,b2]]);print('s={:.2f}\tp={:.2e}'.format(s,p))
        
        p_marker= "{:.1e}".format(p)
        if p_marker[-2]=='0':
            p_marker = p_marker[:-2]+p_marker[-1]
        p_marker='*'
        if p<0.001:
            p_marker = '**'

        if p<0.001:
            plt.text(position,100*b2/b1 , p_marker,ha = 'center',fontsize=14)
        elif p<0.05:
            plt.text(position,100*b2/b1 , p_marker ,ha = 'center',fontsize=14)
        position+=1

    plt.axes().set_xticks(np.arange(position))
    plt.axes().set_xticklabels(xticklabels,rotation=30,ha = 'right',fontsize=16)
#     sns.despine(offset=None, trim=False)
    plt.ylabel(ylabel,fontsize=16)
    plt.xlim([-.5,position-.5])
    plt.ylim(ylim)
#     plt.axes().tick_params(axis='x',direction='out', length=3, width=.8, colors='black')
#     plt.axes().tick_params(axis='y',direction='out', length=3, width=.8, colors='black')
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.02,dpi=600,transparent=True)
    plt.close()



    
def return_test_score(a,b,c):
    p1 = len(a[0])/len(a[1])
    p2 = len(b[0])/len(b[1])
    p3 = len(c[0])/len(c[1])
    print('#gained\t',len(a[0]),len(a[1]),p1)
    print('#lost\t',len(b[0]),len(b[1]),p2)
    print('#const\t',len(c[0]),len(c[1]),p3)
    s,p = stats.fisher_exact([[len(b[0]),len(b[1])-len(b[0])],[len(a[0]),len(a[1])-len(a[0])]])
    print('gained vs. lost:\t',s,p)
    s,p = stats.fisher_exact([[len(c[0]),len(c[1])-len(c[0])],[len(a[0]),len(a[1])-len(a[0])]])
    print('gained vs. const:\t',s,p)



outdir = 'f1_dynamic_gene_enrichment_CUTLL1'
os.makedirs(outdir,exist_ok=True) 
fcthre=np.log2(1.2)
pthre=0.01
flag ='ori'

# intra-domain genes
gained_DomainGenes,lost_DomainGenes,const_DomainGenes = CTCF_TALL_modules.return_intra_domain_genes('T-ALL')
# CUTLL1 vs. CD4
cutll1_all,cutll1_upgenes,cutll1_dngenes = CTCF_TALL_modules.return_CUTLL1_vs_CD4_deg(fcthre=fcthre,pthre=pthre)
#NOTCH1 targets
notch_targets = CTCF_TALL_modules.return_Notch1_target_genes()
# CUTLL1 shCTCF
shCTCF_all,shCTCF_upgenes,shCTCF_dngenes = CTCF_TALL_modules.return_CUTLL1_shCTCF_deg(fcthre=fcthre,pthre=pthre)
# CUTLL1 GSI     
gsi_all_genes, gsi_up_genes, gsi_dn_genes = CTCF_TALL_modules.return_CUTLL1_GSI_deg(fcthre=fcthre,pthre=pthre,flag=flag)
gsi_wo_all_genes, gsi_wo_up_genes, gsi_wo_dn_genes = CTCF_TALL_modules.return_CUTLL1_GSI_wo_deg(fcthre=fcthre,pthre=pthre,flag=flag)       
# dynamic strict vs. tall gained target
gsi_dynamic_genes = np.intersect1d(gsi_dn_genes,gsi_wo_up_genes)
gsi_dynamic_genes = np.union1d(gsi_dn_genes,gsi_wo_up_genes)

#######
## compare of enrichment of dynamic/up genes in gained/lost/const domains
#######

#### % of intra-domain dynamic genes
gained_intra_domain_dynamic = np.intersect1d(gsi_dynamic_genes,gained_DomainGenes)
lost_intra_domain_dynamic = np.intersect1d(gsi_dynamic_genes,lost_DomainGenes)
const_intra_domain_dynamic = np.intersect1d(gsi_dynamic_genes,const_DomainGenes)
a = [gained_intra_domain_dynamic,gained_DomainGenes]
b = [lost_intra_domain_dynamic,lost_DomainGenes]
c = [const_intra_domain_dynamic,const_DomainGenes]
print('\n====\n#intra-domain-dynamic/#intra-domain-genes')
return_test_score(a,b,c)

#### % of intra-domain up-regulated genes
gained_intra_domain_upgenes = np.intersect1d(cutll1_upgenes,gained_DomainGenes)
lost_intra_domain_upgenes = np.intersect1d(cutll1_upgenes,lost_DomainGenes)
const_intra_domain_upgenes = np.intersect1d(cutll1_upgenes,const_DomainGenes)
a = [gained_intra_domain_upgenes,gained_DomainGenes]
b = [lost_intra_domain_upgenes,lost_DomainGenes]
c = [const_intra_domain_upgenes,const_DomainGenes]
print('\n====\n#intra-domain-CUTLL1-up/#intra-domain-genes')
return_test_score(a,b,c)

#### % of intra-domain dynamic up-regulated genes
gained_intra_domain_dynamic_upgenes = np.intersect1d(gained_intra_domain_upgenes,gsi_dynamic_genes)
lost_intra_domain_dynamic_upgenes = np.intersect1d(lost_intra_domain_upgenes,gsi_dynamic_genes)
const_intra_domain_dynamic_upgenes = np.intersect1d(const_intra_domain_upgenes,gsi_dynamic_genes)
a = [gained_intra_domain_dynamic_upgenes,gained_intra_domain_upgenes]
b = [lost_intra_domain_dynamic_upgenes,lost_intra_domain_upgenes]
c = [const_intra_domain_dynamic_upgenes,const_intra_domain_upgenes]
print('\n====\n#intra-domain-dynamic-CUTLL1-up/#intra-domain-CUTLL1-up')
return_test_score(a,b,c)


#######
## compare of enrichment of dynamic/up genes in shCTCF down genes
#######


# GSI
base_lists = [gsi_all_genes,shCTCF_dngenes]
targets = gsi_dn_genes
ylabel='Down-regulated \n genes in GSI (%)'
ylim=[0,6.5]
xticklabels = ['All genes','shCTCF down']
colors = ['silver','firebrick']
figname = outdir+os.sep+'CUTLL1_GSI_DEG_down_enrichIn_shCTCF_down.pdf'
enrichment_compr_bar(base_lists,targets,ylabel,ylim,xticklabels,colors,figname)

# GSI wo
base_lists = [gsi_all_genes,shCTCF_dngenes]
targets = gsi_wo_up_genes
ylabel='Up-regulated\n genes in GSI-w16h (%)'
ylim=[0,12.5]
xticklabels = ['All genes','shCTCF down']
colors = ['silver','firebrick']
figname = outdir+os.sep+'CUTLL1_GSI_w16h_DEG_up_enrichIn_shCTCF_down.pdf'
enrichment_compr_bar(base_lists,targets,ylabel,ylim,xticklabels,colors,figname)

