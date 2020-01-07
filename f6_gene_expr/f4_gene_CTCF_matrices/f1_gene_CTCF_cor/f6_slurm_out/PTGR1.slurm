#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTGR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTH1R 
