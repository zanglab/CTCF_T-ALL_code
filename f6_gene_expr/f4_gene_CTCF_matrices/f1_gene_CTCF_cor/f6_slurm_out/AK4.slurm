#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AK4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AK5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AK6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AK7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AK8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AK9 
