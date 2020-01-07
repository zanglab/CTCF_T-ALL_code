#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HOXB4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXB5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXB6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXB7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXB8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXB9 
