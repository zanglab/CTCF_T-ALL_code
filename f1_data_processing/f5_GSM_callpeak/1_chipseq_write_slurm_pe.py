import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd

#import re,bisect
plus = re.compile('\+')
#minus = re.compile('\-')
name = re.compile('\#name')
treat = re.compile('\#treat')
control = re.compile('\#control')


def write_fq_process(outname,slurmout,genome_index,splitfile1,splitfile2,fq_base,fqdir,outdir,aln=False):
    # write the process of each of the fq file
    fqfile1 = fqdir+os.sep+splitfile1
    fqfile2 = fqdir+os.sep+splitfile2
    saifile1 = outdir+os.sep+outname+'.1.sai'
    saifile2 = outdir+os.sep+outname+'.2.sai'
    samfile = outdir+os.sep+outname+'.sam'
    bamfile = outdir+os.sep+outname+'.bam'
    #bedfile = outdir+os.sep+fq_base+'.bed'
    if aln:
        # reads<=50bp
        slurmout.write('bwa aln  -q 5 -l 32 -k 2 -t 8 {} {} > {}\n'.format(genome_index,fqfile1,saifile1))
        slurmout.write('bwa aln  -q 5 -l 32 -k 2 -t 8 {} {} > {}\n'.format(genome_index,fqfile2,saifile2))
        slurmout.write('bwa sampe {} {} {} {} {} > {}\n'.format(genome_index,saifile1,saifile2,fqfile1,fqfile2,samfile))
    else:
        # reads >=75bp
        slurmout.write('bwa mem -t 8 {} {} {} > {}\n'.format(genome_index,fqfile1,fqfile2,samfile))
    slurmout.write('samtools view -bS -q 1 -@ 8 {} > {}\n\n'.format(samfile,bamfile))
    #slurmout.write('bedtools bamtobed -i {} > {}'.format(bamfile,bedfile))
    return bamfile



def write_slurm(block,infile,slurmdir,fqdir,outdir,species):
    block = [i for i in block.split('\n') if len(i)>0]
    controlline = len(block)
    if len(block)>0:
        #print(block)
        for i in np.arange(len(block)):
            if name.match(block[i]):
                nameline = i
            if treat.match(block[i]):
                treatline = i
            if control.match(block[i]):
                controlline = i
        outname = block[nameline+1]
        treatnames = block[treatline+1:controlline]
        controlnames = block[controlline+1:len(block)]
        outdir = outdir+os.sep+outname
        os.makedirs(outdir,exist_ok=True)

        with open(slurmdir+os.sep+'{}.slurm'.format(outname),'w') as slurmout:
            slurmout.write('''#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
''')

            #######################
            # REMEMBER TO CHANGE THE MEMORY !!!
            #########################
            
            slurmout.write('#SBATCH -o {}/{}.out \n\n'.format(slurmdir,outname))
            slurmout.write('module load bwa\n')
            slurmout.write('module load samtools\n\n')
            
            # set the parameters
            aln=True
            #split_flag='.fastq'
            split_flag1 = '_1.fastq.gz'
            split_flag2 = '_2.fastq.gz'          
            match1 = re.compile(split_flag1)   
            match2 = re.compile(split_flag2)            
            
            # here starts the process of fq files
            genome_index = '/scratch/zw5j/data/index/bwa_index/{}/{}.fa'.format(species,species)
            macs_treat,macs_ctrl = [],[]
            for fqname in treatnames:
                if match1.search(fqname) and not match2.search(fqname):
                    fq_base = match1.split(os.path.basename(fqname))[0] # basename for matched paired end datasets
                    splitfile1 = [i for i in treatnames if (match1.search(i) and re.compile(fq_base).search(i))][0]
                    splitfile2 = [i for i in treatnames if (match2.search(i) and re.compile(fq_base).search(i))][0]
                    bamfile = write_fq_process(outname,slurmout,genome_index,splitfile1,splitfile2,fq_base,fqdir,outdir,aln=aln)
                    macs_treat.append(bamfile)
            for fqname in controlnames:
                if match1.search(fqname) and not match2.search(fqname):
                    fq_base = match1.split(os.path.basename(fqname))[0] # basename for matched paired end datasets
                    splitfile1 = [i for i in controlnames if (match1.search(i) and re.compile(fq_base).search(i))][0]
                    splitfile2 = [i for i in controlnames if (match2.search(i) and re.compile(fq_base).search(i))][0]
                    bamfile = write_fq_process(outname,slurmout,genome_index,splitfile1,splitfile2,fq_base,fqdir,outdir,aln=aln)
                    macs_ctrl.append(bamfile)

            macs_g = 'hs' if re.search('hg',species) else 'mm'
            if len(macs_ctrl)>0:
                slurmout.write('macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g {} -f BAMPE -t \\\n{} -c \\\n{} -n {}/{}\n\n'.format(macs_g,' \\\n'.join(macs_treat),' \\\n'.join(macs_ctrl),outdir,outname))
            else:
                slurmout.write('macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g {} -f BAMPE -t \\\n{} -n {}/{}\n\n'.format(macs_g,' \\\n'.join(macs_treat),outdir,outname))    
            slurmout.write('bdg2bw {}/{}_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/{}/{}.chrom.sizes\n\n'.format(outdir,outname,species,species))
        
            # qc for each block
            qc_out_file = slurmdir+os.sep+'{}_qc.log'.format(outname)
            slurmout.write('python 2_chipseq_qc_pe_macs2_fragments_no_reps.py -i {} -o {} -s {} -b {} -f1 {} -f2 {}'.format(infile,qc_out_file,species,outname,split_flag1,split_flag2)) 




def main(infile,slurmdir,species):
    os.makedirs(slurmdir,exist_ok=True)
    with open(infile) as inf:
        lines = inf.readlines()
        filedir = lines[0].strip().split('\t')[1]
        rundir = lines[1].strip().split('\t')[1]
        names_block = ''
        for line in lines[2:]:                   
            if len(line)>0:               
                if plus.match(line):
                    write_slurm(names_block,infile,slurmdir,filedir,rundir,species)
                    names_block = ''                    
                else:
                    names_block = names_block+line
                    
        write_slurm(names_block,infile,slurmdir,filedir,rundir,species)     
            
   	        




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of match name list', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of slurm files,default: current dir', metavar = '<dir>',default='./')
    parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<3:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile,args.outdir,args.species)
