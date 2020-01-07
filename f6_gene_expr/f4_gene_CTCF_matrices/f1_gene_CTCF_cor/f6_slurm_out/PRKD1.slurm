#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRKD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKDC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKG2 
