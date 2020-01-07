#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ABCA4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCA9 
