#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ITGAD.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGAD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGAE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGAL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGAM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGAV 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGAX 
