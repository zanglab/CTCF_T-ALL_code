#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HEXIM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HEXIM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEXIM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEY2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEYL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HFE 
