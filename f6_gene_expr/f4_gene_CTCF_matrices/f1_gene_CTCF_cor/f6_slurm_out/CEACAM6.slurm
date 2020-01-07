#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CEACAM6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CEACAM6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEACAM7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEACAM8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEBPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEBPB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEBPD 
