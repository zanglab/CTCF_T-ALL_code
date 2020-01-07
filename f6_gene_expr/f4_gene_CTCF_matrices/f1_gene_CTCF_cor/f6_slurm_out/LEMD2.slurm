#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LEMD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LEMD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LEMD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LENEP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LENG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LENG8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LENG9 
