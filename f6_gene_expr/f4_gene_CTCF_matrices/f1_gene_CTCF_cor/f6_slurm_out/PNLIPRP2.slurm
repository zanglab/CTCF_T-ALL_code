#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PNLIPRP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PNLIPRP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNLIPRP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNMA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNMA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNMA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNMA5 
