#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LAPTM5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LAPTM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LARP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LARP1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LARP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LARP4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LARP6 
