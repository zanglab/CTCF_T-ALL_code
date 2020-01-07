#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PUM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PUM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PUM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PURA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PURB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PURG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PUS1 
