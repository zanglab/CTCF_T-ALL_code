#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VAMP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VAMP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VAMP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VAMP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VAMP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VAMP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VAMP7 
