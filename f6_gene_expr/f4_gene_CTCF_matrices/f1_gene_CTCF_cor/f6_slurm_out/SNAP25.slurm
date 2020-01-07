#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SNAP25.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SNAP25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNAP29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNAP47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNAP91 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNAPC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNAPC2 
