#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DUSP6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUX4 
