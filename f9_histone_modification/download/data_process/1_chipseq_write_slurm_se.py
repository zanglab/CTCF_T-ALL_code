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


def write_fq_process(outname,slurmout,genome_index,fqname,fqdir,outdir,aln=False,split_flag='.fastq'):
    # write the process of each of the fq file
    #### the outname is needed to specify the output files
    fq_base = outname+'_'+os.path.basename(fqname).split(split_flag)[0] # remember to add prefix for each block
    fqfile = fqdir+os.sep+fqname
    saifile = outdir+os.sep+fq_base+'.sai'
    samfile = outdir+os.sep+fq_base+'.sam'
    bamfile = outdir+os.sep+fq_base+'.bam'
    #bedfile = outdir+os.sep+fq_base+'.bed'
    if aln:
        # reads<=50bp
        slurmout.write('bwa aln  -q 5 -l 32 -k 2 -t 8 {} {} > {}\n'.format(genome_index,fqfile,saifile))
        slurmout.write('bwa samse {} {} {} > {}\n'.format(genome_index,saifile,fqfile,samfile))
    else:
        # reads >=75bp
        slurmout.write('bwa mem -t 8 {} {} > {}\n'.format(genome_index,fqfile,samfile))
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
#SBATCH -n 8
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
            split_flag='.fastq.gz'
            
            # here starts the process of fq files
            genome_index = '/scratch/zw5j/data/index/bwa_index/{}/{}.fa'.format(species,species)
            macs_treat,macs_ctrl = [],[]
            for fqname in treatnames:
                bamfile = write_fq_process(outname,slurmout,genome_index,fqname,fqdir,outdir,aln=aln,split_flag=split_flag)
                macs_treat.append(bamfile)
            for fqname in controlnames:
                bamfile = write_fq_process(outname,slurmout,genome_index,fqname,fqdir,outdir,aln=aln,split_flag=split_flag)
                macs_ctrl.append(bamfile)
            macs_g = 'hs' if re.search('hg',species) else 'mm'
            if len(macs_ctrl)>0:
                slurmout.write('macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g {} -t \\\n{} -c \\\n{} -n {}/{}\n\n'.format(macs_g,' \\\n'.join(macs_treat),' \\\n'.join(macs_ctrl),outdir,outname))
            else:
                slurmout.write('macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g {} -t \\\n{} -n {}/{}\n\n'.format(macs_g,' \\\n'.join(macs_treat),outdir,outname))    
            slurmout.write('bdg2bw {}/{}_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/{}/{}.chrom.sizes\n\n'.format(outdir,outname,species,species))
        
            # qc for each block
            qc_out_file = slurmdir+os.sep+'{}_qc.log'.format(outname)
            slurmout.write('python 2_chipseq_qc_se.py -i {} -o {} -s {} -b {} -f {}'.format(infile,qc_out_file,species,outname,split_flag)) 




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
