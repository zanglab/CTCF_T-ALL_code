#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AKIRIN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AKIRIN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKNA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKNAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKR1A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKR1B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKR1B10 
