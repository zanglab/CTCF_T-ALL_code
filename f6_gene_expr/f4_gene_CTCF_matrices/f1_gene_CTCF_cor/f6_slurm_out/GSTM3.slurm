#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GSTM3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTP1 
