#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SOWAHA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SOWAHA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOWAHB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOWAHC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOWAHD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX10 
