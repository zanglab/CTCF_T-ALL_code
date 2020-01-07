#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DYNC1I2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DYNC1I2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYNC1LI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYNC1LI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYNC2H1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYNC2LI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYNLL1 
