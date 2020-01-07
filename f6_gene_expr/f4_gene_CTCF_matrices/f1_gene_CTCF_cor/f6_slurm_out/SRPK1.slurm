#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SRPK1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SRPK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRPK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRPK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRPRB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRPX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRPX2 
