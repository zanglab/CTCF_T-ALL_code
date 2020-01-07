#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ASPHD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ASPHD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASPM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASPN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASPRV1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASPSCR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASRGL1 
