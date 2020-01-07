#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM9B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM9B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMIE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMIGD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMIGD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMIGD3 
