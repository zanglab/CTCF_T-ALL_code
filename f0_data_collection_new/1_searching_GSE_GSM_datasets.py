import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
import multiprocessing
#from GenomeData import *
#import scipy
#from scipy import stats
#from scipy.cluster.hierarchy import linkage,dendrogram,cut_tree
#sys.setrecursionlimit(10000)
import urllib.request
from bs4 import BeautifulSoup
import multiprocessing
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
                                

def add_gsm(gse):
    # get info from NCBI for each GSEid
    search_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}".format(gse)#;print(gse)
    search_response = urllib.request.urlopen(search_url)
    search_website = search_response.read()
    search_response.close()
    #search_html = search_website.decode("utf-8")
    soup = BeautifulSoup(search_website,'html.parser')#;print(soup);exit(0)
    lines = soup.find('td').get_text(separator=u'\n').splitlines() #replace <br> with \n
    #print(lines);#exit(0)
    gse_samples = []
    for ele in np.arange(len(lines)):
        # get GSM for CTCF
        if re.match("GSM",lines[ele]):
            print('{}\t{}'.format(lines[ele],gse))
#             for next in np.arange(ele+1,len(lines)): # gsm title listed in next line
#                 if len(lines[next])>0:
#                     gse_samples.append([lines[ele],gse,lines[next]]);break
    return gse_samples


def main(infile):


    # series file download from:
    # https://www.ncbi.nlm.nih.gov/geo/browse/
    # all ChIP-seq, all human
    
    with open('f0_infile/hg_series.csv','rb') as inf1,open('f0_infile/hg_series2.csv','rb') as inf2:
        series1 = pd.read_csv(inf1,index_col=0)
        series2 = pd.read_csv(inf2,index_col=0)
    series = pd.concat([series1,series2])
    del series1,series2

#     series=series.iloc[:10]
#     print(series)
    pool = multiprocessing.Pool(processes=80)
    all_gsm_samples = pool.map_async(add_gsm,series.index,chunksize=1)
    pool.close()
    pool.join()
    
    
    
    
#     gsm_gse = pd.DataFrame()
#     for gsm_samples in all_gsm_samples.get():
#         try:
#             for sample in gsm_samples:
#                 gsm_gse.loc[sample[0],'GSM'] = sample[0]
#                 gsm_gse.loc[sample[0],'GSE'] = sample[1]
# #                 gsm_gse.loc[sample[0],'Title'] = sample[2]
#                 gsm_gse.loc[sample[0],'PMID'] = series.loc[sample[1],'PubMed ID']
#                 gsm_gse.loc[sample[0],'Contact'] = series.loc[sample[1],'Contact']
#                 gsm_gse.loc[sample[0],'ReleaseDate'] = series.loc[sample[1],'Release Date']
#         except:
#             print(gsm_samples,'failed')
#     gsm_gse.to_csv('NCBI_Homo_ChIPseq_GSE_GSM_201906.csv',sep='\t',index=False)



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
