#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANGPTL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANGPTL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANGPTL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANGPTL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANGPTL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANGPTL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANGPTL6 
