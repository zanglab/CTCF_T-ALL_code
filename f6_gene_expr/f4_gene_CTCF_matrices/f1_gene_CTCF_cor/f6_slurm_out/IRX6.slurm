#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IRX6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IRX6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ISCA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ISCA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ISCU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ISG15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ISG20 
