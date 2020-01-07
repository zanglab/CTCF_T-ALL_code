#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AIRE.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AIRE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AJAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AJUBA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AK3 
