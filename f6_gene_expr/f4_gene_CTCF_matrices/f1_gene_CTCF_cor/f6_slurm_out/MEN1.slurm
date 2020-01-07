#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MEN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MEN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEOX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEOX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEP1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEP1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEPCE 
