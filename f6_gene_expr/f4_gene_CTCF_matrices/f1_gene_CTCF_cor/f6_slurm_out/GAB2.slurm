#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GAB2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GAB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABARAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABARAPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABARAPL2 
