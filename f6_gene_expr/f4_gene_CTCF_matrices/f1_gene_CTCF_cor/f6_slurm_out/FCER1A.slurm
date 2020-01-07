#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FCER1A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FCER1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCER1G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCER2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCGBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCGR1A 
