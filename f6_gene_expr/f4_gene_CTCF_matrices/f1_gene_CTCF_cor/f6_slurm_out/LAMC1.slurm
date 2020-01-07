#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LAMC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMP3 
