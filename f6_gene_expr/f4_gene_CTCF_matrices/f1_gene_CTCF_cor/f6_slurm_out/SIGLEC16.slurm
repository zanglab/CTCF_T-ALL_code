#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SIGLEC16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SIGLEC16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIGLEC17P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIGLEC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIGLEC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIGLEC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIGLEC8 
