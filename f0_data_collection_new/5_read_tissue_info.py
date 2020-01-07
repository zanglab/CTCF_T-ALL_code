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
                                

def read_cellline(gsm,description,cellline_tissue_df):
    
    mapped_celllines, mapped_tissues, cancerous, treatment = set(),set(),set(),set()
    celltype = ''
    #print(description)
    description = re.split(':|\s+',description)
    description = [''.join(re.split('_|-|\.',ele)) for ele in description]
    #print(description);exit(0)
    
    for ele in description:
        for match_cellline in cellline_tissue_df.index:
            match_cellline_character = ''.join(re.split('_|-|\.',match_cellline))
            if re.search(match_cellline_character,ele,flags=re.I):
                mapped_celllines.add(match_cellline)
                mapped_tissues.add(cellline_tissue_df.loc[match_cellline,'tissue'])
                if cellline_tissue_df.loc[match_cellline,'cancerous']!='NA':
                    cancerous.add(cellline_tissue_df.loc[match_cellline,'cancerous'])
                if cellline_tissue_df.loc[match_cellline,'celltype']!='NA':
                    celltype= cellline_tissue_df.loc[match_cellline,'celltype']

                
        for tissue in ['CD4','breast','colon','kidney','liver','lung','pancreas','prostate']:
            if re.search(tissue,ele,flags=re.I):
                mapped_tissues.add(tissue)
        
        for cancerous_type in ['cancer','cancerous','carcinoma']:
            if re.search(cancerous_type,ele,flags=re.I):
                cancerous.add('Y')

               
        for control_type in ['control','input','IgG']:
            if re.search(control_type,ele,flags=re.I):
                treatment.add(control_type)    
        if re.search('CTCFL',ele):
            print(ele,gsm)    
    mapped_celllines, mapped_tissues,cancerous, treatment = ', '.join(mapped_celllines), ', '.join(mapped_tissues),', '.join(cancerous),', '.join(treatment)
    #mapped_celllines, mapped_tissues,cancerous, treatment = ', '.join(list(mapped_celllines)[-1:]), ', '.join(mapped_tissues),', '.join(cancerous),', '.join(treatment)
    
    return mapped_celllines,mapped_tissues,celltype,cancerous,treatment


def main(infile):
    # this is to check if the 377 CTCF datasets in the newly collected datasets
    
    ctcf_gsm = pd.read_csv('NCBI_Homo_ChIPseq_GSE_GSM_CTCF_filtered.csv',sep='\t')
    ctcf_gsm.columns = ['GSM', 'GSE', 'antibody', 'characteristics', 'library','organism', 'sourcename', 'title']
    # add release data/PubMed ID
    with open('f0_infile/hg_series.csv','rb') as inf1,open('f0_infile/hg_series2.csv','rb') as inf2:
        series1 = pd.read_csv(inf1,index_col=0)
        series2 = pd.read_csv(inf2,index_col=0)
    series = pd.concat([series1,series2])
    series = series[['PubMed ID', 'Contact', 'Release Date']]
    series.columns = ['PMID', 'Contact', 'ReleaseDate']
#     print(series)
    ctcf_gsm = ctcf_gsm.merge(series,how='left',left_on=["GSE"],right_on=series.index)
    
#     exit()
    cellline_tissue_df = pd.read_excel('f0_infile/cellline_tissue_201906.xlsx',index_col=0)
    cellline_tissue_df = cellline_tissue_df.fillna('NA')
    # get the tissue info for each label
    for gsm_index in ctcf_gsm.index:
        title = ctcf_gsm.loc[gsm_index,'title']
        sourcename = ctcf_gsm.loc[gsm_index,'sourcename']
        characteristics = ctcf_gsm.loc[gsm_index,'characteristics']
        description = ' && '.join([title,sourcename,characteristics])
        cellline,tissue,celltype,cancerous,treatment = read_cellline(ctcf_gsm.loc[gsm_index,'GSM'],description,cellline_tissue_df)
        ctcf_gsm.loc[gsm_index,'cellline'] = cellline
        ctcf_gsm.loc[gsm_index,'tissue'] = tissue
        ctcf_gsm.loc[gsm_index,'celltype'] = celltype
        ctcf_gsm.loc[gsm_index,'cancerous'] = cancerous
        ctcf_gsm.loc[gsm_index,'treatment'] = treatment
    
    ctcf_gsm = ctcf_gsm[['GSM','GSE','cellline','tissue','celltype','cancerous','treatment',\
        'PMID', 'Contact', 'ReleaseDate', 'antibody','library', 'organism',\
        'title', 'sourcename', 'characteristics']]
    ctcf_gsm.to_csv('CTCF_ChIPseq_GSE_GSM_description_addTissue.csv',sep='\t',index=False)


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
