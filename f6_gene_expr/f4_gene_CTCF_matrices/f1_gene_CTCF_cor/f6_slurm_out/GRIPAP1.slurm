#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GRIPAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIPAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRK5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRK6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRK7 
