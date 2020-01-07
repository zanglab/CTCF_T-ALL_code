#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYOD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYOD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYOF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYOG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYOM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYOM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYOM3 
