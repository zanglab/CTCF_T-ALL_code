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
                                

def gsm_get_info(gsm):
    #print(ctcf_df.loc[label]);exit(0)
    antibody_keywords = "antibody"
    title_keywords = "Title"
    source_keywords = "Source name"
    organism_keywords = "Organism"
    characteristics_keywords = "Characteristics"
    library_keywords = "Library strategy"
    antibody,title,sourcename,organism,characteristics,library = [],[],[],[],[],[] #init as none for each id  
    
    # get info from NCBI for each GSMid
    gsm_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}".format(gsm);print(gsm)
    gsm_response = urllib.request.urlopen(gsm_url)
    gsm_website = gsm_response.read()
    gsm_response.close()
    #gsm_html = gsm_website.decode("utf-8")
    soup = BeautifulSoup(gsm_website,'html.parser')#;print(soup);exit(0)
    lines = soup.find('td').get_text(separator=u'\n').splitlines() #replace <br> with \n
    #print(lines);exit(0)
    for ele in np.arange(len(lines)):
        # Title 
        if re.match(title_keywords,lines[ele]):
            for next in np.arange(ele+1,len(lines)): # source name listed in next line
                if len(lines[next])>0:
                    title.append(lines[next]);break

        # Source name
        if re.match(source_keywords,lines[ele]):
            for next in np.arange(ele+1,len(lines)): # source name listed in next line
                if len(lines[next])>0:
                    sourcename.append(lines[next]);break

        # organism
        if re.match(organism_keywords,lines[ele]):
            for next in np.arange(ele+1,len(lines)): # source name listed in next line
                if len(lines[next])>0:
                    organism.append(lines[next]);break

        # library
        if re.match(library_keywords,lines[ele]):
            for next in np.arange(ele+1,len(lines)): # source name listed in next line
                if len(lines[next])>0:
                    library.append(lines[next]);break

        # characteristics
        if re.match(characteristics_keywords,lines[ele]):
            for next in np.arange(ele+1,len(lines)): # source name listed in next line
                if len(lines[next])>0:
                    characteristics.append(lines[next]);#print(lines[next])
                    if len(lines[next+1])==0:
                        break 
    # for each term in characteristics, check if exist antibody info
    for line in characteristics:
        term = line.split(':');#print(term)
        if re.search('antibody',term[0],flags=re.I):
            antibody.append(term[1])
    
    title = ' && '.join(title)
    characteristics = ' && '.join(characteristics) 
    sourcename = ' && '.join(sourcename) 
    organism = ' && '.join(organism)
    library = ' && '.join(library)
    antibody = ' && '.join(antibody)
    #print(title,sourcename,organism,characteristics,library,antibody);exit(0)               
    with open('f2_gsmInfo/{}.txt'.format(gsm),'w') as outf:
        outf.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gsm,antibody,title,sourcename,organism,characteristics,library))
    return gsm,antibody,title,sourcename,organism,characteristics,library


def main(infile):

    ctcf_gsm = pd.read_csv('NCBI_Homo_ChIPseq_GSE_GSM_201906.csv',sep='\t',index_col=0,header=None)
#     print(ctcf_gsm);exit()
    exist_data = pd.read_csv('NCBI_Homo_ChIPseq_GSE_GSM_allInfo.csv',sep='\t',index_col=0)
    ctcf_gsm = ctcf_gsm.loc[ctcf_gsm.index.difference(exist_data.index)]
#     ctcf_gsm = ctcf_gsm.iloc[:2000]   
    
    # get the tissue info for each label
    pool = multiprocessing.Pool(processes=10)
    all_gsm_descriptions = pool.map_async(gsm_get_info,ctcf_gsm.index,chunksize=1)
    pool.close()
    pool.join()
    






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
