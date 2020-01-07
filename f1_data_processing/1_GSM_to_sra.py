import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
# import multiprocessing
#from GenomeData import *
#import scipy
#from scipy import stats
#from scipy.cluster.hierarchy import linkage,dendrogram,cut_tree
#sys.setrecursionlimit(10000)
# import urllib.request
# from bs4 import BeautifulSoup
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
import urllib.request
from bs4 import BeautifulSoup                                

def GSM_to_SRR(GSM_id):
    # fetch the SRR from for GSM id
    gsm_keywords = "sra\?term"
    gsm_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}".format(GSM_id)
    gsm_response = urllib.request.urlopen(gsm_url)
    gsm_website = gsm_response.read()
    gsm_response.close()
    gsm_html = gsm_website.decode("utf-8")
    #grep the SRX id
    for i in gsm_html.splitlines():
        if re.search(gsm_keywords,i):
            SRX_lines=re.split('=|"',i)
            for SRX_line in SRX_lines:
                if SRX_line.startswith('SRX'):
                   SRX_id = SRX_line
#     print(SRX_id)                            
    srx_keywords = "run=SRR"
    srx_url = "https://www.ncbi.nlm.nih.gov/sra?term={}".format(SRX_id)
    srx_response = urllib.request.urlopen(srx_url)
    srx_website = srx_response.read()
    srx_response.close()
    srx_html = srx_website.decode("utf-8")
    #grep all the SRR ids
    for i in srx_html.splitlines():
        #print(i)
        if re.search(srx_keywords,i):
            SRR_lines=re.split('run=|"',i)
            SRR_ids = [SRR_line for SRR_line in SRR_lines if SRR_line.startswith('SRR')]  
#     print(SRR_ids)
    gsm_url = "https://www.ncbi.nlm.nih.gov/sra?term={}".format(SRX_id)
    gsm_response = urllib.request.urlopen(gsm_url)
    gsm_website = gsm_response.read()
    gsm_response.close()
    #gsm_html = gsm_website.decode("utf-8")
    soup = BeautifulSoup(gsm_website,'html.parser')#;print(soup);exit(0)
    lines = soup.find('div').get_text(separator=u'\n').splitlines() #replace <br> with \n
    for ele in np.arange(len(lines)):
        # Title 
        if re.match("Layout",lines[ele]):
            layout = lines[ele+1];break
#     print(layout)
    return SRR_ids,layout



def main(infile):

    outdir='f1_GSM2SRR'
    os.makedirs(outdir,exist_ok=True)
     
    df = pd.read_csv('/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f0_data_collection_new/hg_occupancy_ChIP/CTCF_ChIPseq_GSE_GSM_description_addTissue_withhJournal_peaks.csv',sep=',',index_col=0)
    df = df[df['peak_nums'].isnull()==True]
    print('#new data:\t',df.shape)
    df_out = pd.DataFrame(columns=['SRR','layout'])
    for gsmid in df.index:
#         print(gsmid)
        try:
            srrs,layout = GSM_to_SRR(gsmid)  
            df_out.loc[gsmid,'SRR'] = '&&'.join(srrs)
            df_out.loc[gsmid,'layout'] = layout
        except:
            print(gsmid,'Not accessible')
    df_out.to_csv(outdir+os.sep+'ctcf_new_GSM_to_SRR.csv')
    #exit() 






    
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
