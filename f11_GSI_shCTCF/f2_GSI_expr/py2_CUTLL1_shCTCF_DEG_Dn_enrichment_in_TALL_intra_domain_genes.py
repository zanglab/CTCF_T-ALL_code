import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.6)
sns.set_style("whitegrid", {'axes.grid' : False})
from scipy import stats
#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
import CTCF_TALL_modules_new

def bar_plot_enrichment(vals,xticklabels,figname):
    a1,a2,b1,b2,c1,c2 = vals
    plt.figure(figsize=(2.2,3.2))
    plot_positions = [1,2,3]
    width = 0.6
    plot_vals = [100*a2/a1,100*b2/b1,100*c2/c1]
    p=plt.bar(plot_positions,plot_vals,width=width,color = ['silver','k','firebrick'])
    ###
    s,p = stats.fisher_exact([[a1-a2,a2],[b1-b2,b2]]);print('Up:',s,p)
    p_label='{:.1e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    plt.text(plot_positions[1]+.6 , plot_vals[1]*1.02 , p_label ,ha = 'right',fontsize=16)
    
    ###
    s,p = stats.fisher_exact([[a1-a2,a2],[c1-c2,c2]]);print('CTCF Up:',s,p)
    p_label='{:.1e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    plt.text(plot_positions[2]+.6 , plot_vals[2]*1.02 , p_label ,ha = 'right',fontsize=16)
    plt.axes().set_xticks(plot_positions)
    plt.axes().set_xticklabels(xticklabels,rotation=25,ha = 'right',fontsize=16)
    sns.despine(offset=None, trim=False)
    plt.ylabel('% of genes down \nin CUTLL1 shCTCF',fontsize=18)
    plt.xlim([0.3,3.5])
#     plt.ylim([0,12])
    plt.axes().tick_params(axis='x',direction='out', length=3, width=.8, colors='black')
    plt.axes().tick_params(axis='y',direction='out', length=3, width=.8, colors='black')
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.02,dpi=600,transparent=True)
    plt.close()



def return_logFC_padj(deseq2_file,fcthre,pthre,flag):
    with open(deseq2_file) as inf:
        df = pd.read_csv(inf,index_col=0)
    upgenes = df[(df['log2FoldChange']>fcthre)&(df['padj']<pthre)].index.values
    dngenes = df[(df['log2FoldChange']<-1*fcthre)&(df['padj']<pthre)].index.values
    print(flag,fcthre,pthre,'up:',len(upgenes),'dn:',len(dngenes))
    return df.index,upgenes,dngenes

def return_ctcf_intra_domain_genes(binding_type,colname):
    # intra domain genes of CTCF bindings
    ctcf_domain_gene_dic = CTCF_TALL_modules_new.return_union_domain_gene()
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    
    tall_gained_intra_domain_genes = np.array([])
    for gained_id in gained_df.index:
        domain_genes = ctcf_domain_gene_dic[str(gained_id)]['intra_domain_genes']
        tall_gained_intra_domain_genes = np.append(tall_gained_intra_domain_genes,domain_genes)
    return tall_gained_intra_domain_genes
    
    

def main(infile):

    outdir='f2_figs'
    os.makedirs(outdir,exist_ok=True)
    
    # deseq2 info
    fcthre,pthre=np.log2(1.2),0.001
    deseq2_CUTLL1_CD4 = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/f4_CTCF_cor_gene_cancer_vs_normal/f2_salmon_deseq2_pca/f3_deseq_out/treated_CUTLL1_vs_ctrl_CD4.csv"
    cutll1_all,cutll1_upgenes,cutll1_dngenes = return_logFC_padj(deseq2_CUTLL1_CD4,fcthre,pthre,'CUTLL1_CD4')
    deseq2_CUTLL1_shCTCF = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/18_MYC-ChIP_shCTCF-RNA/shCTCF_RNA/f0_processing/salmon_Deseq2/salmon_Deseq2_pca/CUTLL1/f3_deseq_out/treated_CUTLL1_shCTCF_vs_ctrl_CUTLL1_PIG.csv"
    shCTCF_all,shCTCF_upgenes,shCTCF_dngenes = return_logFC_padj(deseq2_CUTLL1_shCTCF,fcthre,pthre,'shCTCF')

    # get the intra domain genes of gained/lost CTCF bindings
    tall_gained_intra_domain_genes = return_ctcf_intra_domain_genes('gained','intra_domain_genes')
    tall_gained_intra_domain_shCTCF_down_genes = set(tall_gained_intra_domain_genes).intersection(shCTCF_dngenes)
    tall_gained_intra_domain_up_genes = set(tall_gained_intra_domain_genes).intersection(cutll1_upgenes)
    tall_gained_intra_domain_up_shCTCF_DOWN_genes = set(tall_gained_intra_domain_up_genes).intersection(shCTCF_dngenes)
    cutll1_upgene_shCTCF_DOWN_genes = set(cutll1_upgenes).intersection(shCTCF_dngenes)
    
    a1 = len(set(shCTCF_all).intersection(cutll1_all))
    a2 = len(set(shCTCF_dngenes).intersection(cutll1_all))
    print('\nshCTCF_all',a1)
    print('shCTCF_dngenes',a2)
    
    b1 = len(cutll1_upgenes)
    b2 = len(cutll1_upgene_shCTCF_DOWN_genes)
    print('\ncutll1_upgenes',b1)
    print('cutll1_upgene_shCTCF_DOWN_genes',b2)

    d1 = len(tall_gained_intra_domain_genes)
    d2 = len(tall_gained_intra_domain_shCTCF_down_genes)
    print('\ntall_gained_intra_domain_genes',d1)
    print('tall_gained_intra_domain_shCTCF_down_genes',d2)# plot logFC
    print(sorted(tall_gained_intra_domain_shCTCF_down_genes))

    c1 = len(tall_gained_intra_domain_up_genes)
    c2 = len(tall_gained_intra_domain_up_shCTCF_DOWN_genes)
    print('\ntall_gained_intra_domain_up_genes',c1)
    print('tall_gained_intra_domain_up_shCTCF_DOWN_genes',c2)# plot logFC
    print(sorted(tall_gained_intra_domain_up_shCTCF_DOWN_genes))

    figname = '{}/CUTLL1_shCTCF_Dn_enrichment_bar.png'.format(outdir)
    xticklabels=['All','T-ALL gained intra-domain','T-ALL gained targets']
    vals = (a1,a2,d1,d2,c1,c2)
    bar_plot_enrichment(vals,xticklabels,figname)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
    parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile)
