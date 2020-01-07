#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GRIK4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIK5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIN2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIN2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIN2C 
