#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AIG1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AIG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AIM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AIMP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AIMP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AIPL1 
