#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o REEP6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i REEP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REG1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REG1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REG3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REG3G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REG4 
