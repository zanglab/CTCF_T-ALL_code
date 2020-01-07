#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HES6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HES6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HES7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HESX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEXA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEXB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEXDC 
