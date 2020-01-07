#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IBA57.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IBA57 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IBSP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IBTK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ICA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ICA1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ICAM1 
