#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SOCS5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SOCS5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOCS6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOCS7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOD3 
