#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o XIAP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i XIAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XIRP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XIRP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XKR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XKR4 
