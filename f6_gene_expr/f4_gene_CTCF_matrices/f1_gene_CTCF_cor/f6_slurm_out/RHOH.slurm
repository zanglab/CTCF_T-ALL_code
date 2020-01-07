#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RHOH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOJ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOQ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOU 
