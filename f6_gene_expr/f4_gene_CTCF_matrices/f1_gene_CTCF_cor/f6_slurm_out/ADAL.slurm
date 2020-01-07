#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADAL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM17 
