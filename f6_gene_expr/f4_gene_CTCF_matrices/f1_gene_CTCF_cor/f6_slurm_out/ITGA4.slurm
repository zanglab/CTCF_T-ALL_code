#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ITGA4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGA9 
