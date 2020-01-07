import os,sys,argparse,glob,re
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import gridspec
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=12
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.6)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
import CTCF_TALL_modules_new
import scipy
from scipy.interpolate import interpn

def scatter_plot(ctcf_df,atac_df,figname):

    fig = plt.figure(figsize=(3.5,3.5))

    x,y = ctcf_df[1].values,atac_df[1].values 
    data, x_e,y_e = np.histogram2d(x,y,bins=20)
    z = interpn(( 0.5*(x_e[1:]+x_e[:-1]),0.5*(y_e[1:]+y_e[:-1]) ),data,np.vstack([x,y]).T, method = 'splinef2d', bounds_error=False )
    idx = z.argsort()
    x,y,z = x[idx],y[idx],z[idx]
    g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=1,marker='o')
    plt.xlabel('CTCF log10(RPKM+0.1)',fontsize=15)
    plt.ylabel('ATAC log10(RPKM+0.1',fontsize=15)
    plt.xlim([-1,2])
    plt.ylim([-1,2])
    plt.savefig(figname+'__scatter.png',bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.close()


def hist_plot(ctcf_df,atac_df,figname):

    x,y = ctcf_df[1].values,atac_df[1].values

    fig = plt.figure(figsize=(3.5,3.5))
    sns.distplot(x)
    plt.title('CTCF')
    plt.xlim([-1,2])
    plt.savefig(figname+'__hist_x.png',bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.close()
    
    fig = plt.figure(figsize=(3.5,3.5))
    sns.distplot(y)
    plt.title('ATAC')
    plt.xlim([-1,2])
    plt.savefig(figname+'__hist_y.png',bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.close()

    
def compr_binding_pattern(ctcf_id_file,ctcf_file,atac_file,binding_pattern_type,prename,outdir):

    with open(ctcf_file) as ctcf_inf, open(atac_file) as atac_inf:
        ctcf_df = pd.read_csv(ctcf_inf,sep='\t',index_col=0,header=None)
        atac_df = pd.read_csv(atac_inf,sep='\t',index_col=0,header=None)
    ctcf_df = np.log10(ctcf_df+.1)
    atac_df = np.log10(atac_df+.1)
    
    # read the CTCF ids for scatter plot
    with open(ctcf_id_file) as ctcf_id_inf:
        ctcf_id = pd.read_csv(ctcf_id_inf,sep='\t',index_col=3,header=None)
    
    # use 1e-6 as cutoff to select CTCF motif regions
    ctcf_id = ctcf_id[ctcf_id[6]<1e-6]
#     print(ctcf_id);print(binding_pattern_type, prename, ctcf_id.shape);exit()
    
    ctcf_df = ctcf_df.loc[ctcf_id.index]
    atac_df = atac_df.loc[ctcf_id.index]

    figname = outdir+os.sep+'ctcf_cor_atac_{}_{}'.format(binding_pattern_type,prename,)
    ctcf_id.to_csv(figname+'.csv',header=None,sep='\t')
    scatter_plot(ctcf_df,atac_df,figname)
    hist_plot(ctcf_df,atac_df,figname)
    
        
    
    
def main():        
            

    outdir = 'f2_scatter_by_Jurkat_CTCF_by_motif'
    os.makedirs(outdir,exist_ok=True)

    data_dir='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f5_TCGA_ATAC/f7_ATAC_CTCF_cor/f0_union_binding_RPKM_pattern/data'
    binding_dir='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f5_TCGA_ATAC/f7_ATAC_CTCF_cor/f0_union_binding_RPKM_pattern/binding_csv'

    sub_union= data_dir+os.sep+'union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped.bed'
    sub_union_JurkatCTCF= data_dir+os.sep+'union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped_Jurkat_CTCF_overlapped.bed'
    sub_union_noJurkatCTCF= data_dir+os.sep+'union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped_Jurkat_CTCF_NOT_overlapped.bed'
    
    
    for binding_pattern_type in ['e200','union_site']:
        ctcf_file = binding_dir+os.sep+'ctcf_RPKM_{}.csv'.format(binding_pattern_type)
        atac_file = binding_dir+os.sep+'atac_RPKM_{}.csv'.format(binding_pattern_type)
        compr_binding_pattern(sub_union,ctcf_file,atac_file,binding_pattern_type,'NoDNAme_withMotif',outdir)
        compr_binding_pattern(sub_union_JurkatCTCF,ctcf_file,atac_file,binding_pattern_type,'NoDNAme_withMotif_JurkatCTCF',outdir)
        compr_binding_pattern(sub_union_noJurkatCTCF,ctcf_file,atac_file,binding_pattern_type,'NoDNAme_withMotif_NoJurkatCTCF',outdir)









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
  
    main()

