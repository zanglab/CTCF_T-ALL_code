#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GNAT3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GNAT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNAZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNB1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNB3 
