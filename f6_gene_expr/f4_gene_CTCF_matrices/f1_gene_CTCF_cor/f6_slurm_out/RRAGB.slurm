#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RRAGB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RRAGB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRAGC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRAGD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRAS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRAS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRBP1 
