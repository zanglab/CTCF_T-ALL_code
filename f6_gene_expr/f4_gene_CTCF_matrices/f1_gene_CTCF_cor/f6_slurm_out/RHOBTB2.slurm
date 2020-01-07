#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RHOBTB2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOBTB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOBTB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOG 
