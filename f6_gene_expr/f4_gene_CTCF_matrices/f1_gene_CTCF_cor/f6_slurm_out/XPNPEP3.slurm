#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o XPNPEP3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i XPNPEP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XPO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XPO4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XPO5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XPO6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XPO7 
