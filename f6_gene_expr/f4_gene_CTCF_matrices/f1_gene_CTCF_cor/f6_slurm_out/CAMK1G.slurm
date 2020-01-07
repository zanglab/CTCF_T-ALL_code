#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CAMK1G.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMK1G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMK2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMK2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMK2D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMK2G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMK2N1 
