#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ETFA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ETFA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ETFB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ETFDH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ETHE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ETNK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ETNK2 
