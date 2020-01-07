#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ASB11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB16 
