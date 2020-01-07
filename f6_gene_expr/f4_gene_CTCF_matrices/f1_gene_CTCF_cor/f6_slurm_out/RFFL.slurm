#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RFFL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RFFL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFNG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFPL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFPL3 
