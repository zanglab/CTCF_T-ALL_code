#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHMP7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHMP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHODL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHORDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHP1 
