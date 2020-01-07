#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LIPJ.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LIPJ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIPK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIPM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIPN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIPT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIPT2 
