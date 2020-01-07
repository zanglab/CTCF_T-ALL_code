#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LCA5L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LCA5L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCE1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCE1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCE1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCE1D 
