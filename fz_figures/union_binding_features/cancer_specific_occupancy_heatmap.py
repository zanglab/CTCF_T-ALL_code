import os,sys,argparse,glob,re
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import gridspec
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
# matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
# matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
import CTCF_TALL_modules_new
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"



def plot_num_data(df,columns,cancertype,indir,outdir,bins=20,xrange=False,yrange=False,logscale=True,base=0.1):
    
    x_col = df[columns[0]]
    y_col = df[columns[1]]#;print(y_col);exit()
    xset = [int(i) for i in np.arange(xrange[0],xrange[1])]
    yset = [int(i) for i in np.arange(yrange[0],yrange[1])]
    
    figname = indir+os.sep+'{}_vs_others_occurrence'.format(cancertype)
    if cancertype=="PRAD_TissueAdded":
        cancertype="PRAD"
    xlabel='Occupancy score in\n {} {} ChIP-seq datasets'.format(xrange[1]-1,cancertype)
    ylabel='Occupancy frequency in\n other {} ChIP-seq datasets'.format(yrange[1]-1)
    title = cancertype
    
#     plot_df = pd.DataFrame(index = np.arange(1,bins+1),columns = xset)
#     plot_df[:] = 0
#     for index in df.index:
#         #x,y = df.loc[index,columns[0]],df.loc[index,columns[1]]
#         x,y = df.loc[index,columns[0]],df.loc[index,columns[1]]-df.loc[index,columns[0]]
#         y_position = int(max(np.ceil(bins*(y-min(yset))/(max(yset)-min(yset))),1))
#         plot_df.loc[y_position,x] = plot_df.loc[y_position,x]+1
#     plot_df = plot_df.fillna(0)
#     plot_df.to_csv(figname+'.csv')
    
    plot_df = pd.read_csv(figname+'.csv',index_col=0)
    plt.figure(figsize=(4.1,3.3))

    if logscale:
        dfmax = max(np.log10(plot_df+base).max())
        base = 0.1**((dfmax*3)/4)
        #print(dfmax,base)
        g = sns.heatmap(np.log10(np.transpose(plot_df)+base),cmap = plt.cm.GnBu,cbar_kws={"shrink":.62})
        #for y in range(plot_df.shape[0]):
        #    for x in range(plot_df.shape[1]):
        #        if 0<plot_df.iloc[plot_df.shape[0]-1-y, x]:
        #            plt.text(x + 0.5 , y + 0.5, '%.0f' % plot_df.iloc[plot_df.shape[0]-1-y, x], #data[y,x] +0.05 , data[y,x] + 0.05
        #            #plt.text(x + 0.5 , y + 0.5, '14444', #data[y,x] +0.05 , data[y,x] + 0.05
        #            horizontalalignment='center',
        #            verticalalignment='center',
        #            fontsize=6,color='w')

        cbar = g.collections[0].colorbar        
        cbar.set_clim(-1*dfmax/2,dfmax)
        #base = 0.1**((dfmax*0.75)
        cbar.set_ticks([np.log10(base)]+[i for i in np.arange(1,dfmax,2)])
        cbar.set_ticklabels(['$0$']+['$10^{}$'.format(int(i)) for i in np.arange(1,dfmax,2)])
    else:
        g = sns.heatmap(plot_df,yticklabels = y_tickets,annot=True,fmt="d", cmap = plt.cm.GnBu,cbar_kws={"shrink":.62})

    plt.axes().set_title('CTCF sites distribution',fontsize=16)
    plt.axes().set_xlabel(ylabel,fontsize=16)
    plt.axes().set_ylabel(xlabel,fontsize=16)
    
    #y_tickets = ['{}%'.format(int(100*(i+1)/bins)) for i in np.arange(bins)]
    plt.xticks([i for i in np.arange(0,bins+1,5)],['{}%'.format(int(100*(i)/bins)) for i in np.arange(0,bins+1,5)],rotation=0) 
    if len(xset)>15:
#         plt.yticks([i+0.5 for i in xset if (len(xset)-1-i)%5==0],[len(xset)-1-i for i in xset if (len(xset)-1-i)%5==0],rotation=0)
        plt.yticks([i+0.5 for i in xset if (i)%5==0],[i for i in xset if (i)%5==0],rotation=0)
        #plt.gca().invert_yaxis()
    else:
        plt.yticks(rotation=0)
    plt.gca().invert_yaxis()   
    plt.savefig(outdir+os.sep+os.path.basename(figname)+'.pdf',bbox_inches = 'tight',pad_inches = 0.1,dpi=600,transparent=True)
    plt.close()
    
    
def main():

    
    indir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f1_cancerType_occupancy/f3_cancer_normal_occurrence_csv/'
    outdir = 'figures/cancer_specific_occupancy_heatmap'
    os.makedirs(outdir,exist_ok=True)
    
    co_occurrence_dir='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f1_cancerType_occupancy/f1_cancerType_binding_occurrence/'
    cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']
    
    for cancertype in cancertypes:
        co_occurrence_file= co_occurrence_dir+os.sep+'{}_binding_occurrence.csv'.format(cancertype)
        with open(co_occurrence_file) as co_occurrence_inf:
            df = pd.read_csv(co_occurrence_file,index_col=0)
            #help(CTCF_TALL_modules);exit()
            kept_ids = CTCF_TALL_modules_new.return_occupancy_filtered()
            df = df.loc[kept_ids.index]
            
            columns = ['cancer_peaks','all_peaks']
            num_cancer_datasets = df.iloc[0,:]['cancer_total']
            xrange=[0,num_cancer_datasets+1]
            yrange=[0,771-num_cancer_datasets+1]
            plot_num_data(df,columns,cancertype,indir,outdir,bins=10,xrange=xrange,yrange=yrange)

#             exit()
            
            
            
            
        

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

