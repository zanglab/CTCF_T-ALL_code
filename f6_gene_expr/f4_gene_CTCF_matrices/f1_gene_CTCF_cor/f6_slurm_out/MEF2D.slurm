#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MEF2D.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MEF2D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEFV 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEGF10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEGF11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEGF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEGF8 
