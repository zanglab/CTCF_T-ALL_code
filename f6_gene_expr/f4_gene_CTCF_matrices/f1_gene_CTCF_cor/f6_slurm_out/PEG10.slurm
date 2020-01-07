#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PEG10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PEG10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PELI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PELI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PELI3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PELO 
