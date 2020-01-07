#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DOCK9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DOCK9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOHH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOK4 
