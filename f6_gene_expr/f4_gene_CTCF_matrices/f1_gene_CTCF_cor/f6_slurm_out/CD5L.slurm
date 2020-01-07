#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CD5L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CD5L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD63 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD68 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD69 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD7 
