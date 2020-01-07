#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TCEA2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TCEA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCEA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCEAL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCEAL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCEAL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCEAL4 
