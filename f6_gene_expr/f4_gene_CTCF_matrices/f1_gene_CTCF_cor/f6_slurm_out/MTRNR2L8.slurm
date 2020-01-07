#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MTRNR2L8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MTRNR2L8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTRNR2L9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTRR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTSS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTSS1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTTP 
