#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANGPTL7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANGPTL7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANHX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKAR 
