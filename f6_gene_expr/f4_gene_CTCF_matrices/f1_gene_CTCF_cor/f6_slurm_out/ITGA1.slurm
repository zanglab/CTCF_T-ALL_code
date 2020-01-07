#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ITGA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGA10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGA11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGA2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGA3 
