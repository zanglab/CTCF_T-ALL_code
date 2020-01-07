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
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
                                

def pmid_get_if(pmid):
    try:
        pmid = pmid.split(';')[0]
    except:
        pass
    # get info from NCBI for each GSMid
    journal=''
    gsm_url = "https://www.ncbi.nlm.nih.gov/pubmed/{}".format(pmid)
    gsm_response = urllib.request.urlopen(gsm_url)
    gsm_website = gsm_response.read()
    gsm_response.close()
    gsm_html = gsm_website.decode("utf-8")
    for line in gsm_html.splitlines():
        if re.search('journal=',line):
            search_lines=re.split('"',line)
            for ele in np.arange(len(search_lines)):
                if re.search('journal=',search_lines[ele]):
                    journal = (search_lines[ele+1]);break
    print(journal)
    return journal


def main(infile):

#     ctcf_gsm = pd.read_excel('CTCF_NCBI_ALLCTCF_manu_revised_20180405.xlsx',sep='\t',index_col=0)
    ctcf_gsm = pd.read_csv('CTCF_ChIPseq_GSE_GSM_description_addTissue.csv',sep='\t')
    for index in ctcf_gsm.index:
        if not pd.isnull(ctcf_gsm.loc[index,'PMID']):
            print(ctcf_gsm.loc[index,'PMID'])
            journal = pmid_get_if(ctcf_gsm.loc[index,'PMID'])
            ctcf_gsm.loc[index,'journal']=journal
    ctcf_gsm.to_csv('CTCF_ChIPseq_GSE_GSM_description_addTissue_withhJournal.csv',sep='\t',index=False)
    
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
