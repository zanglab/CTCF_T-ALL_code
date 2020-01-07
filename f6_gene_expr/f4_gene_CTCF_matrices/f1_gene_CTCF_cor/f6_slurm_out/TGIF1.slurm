#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TGIF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TGIF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGIF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGIF2-C20orf24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGIF2LX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGM2 
