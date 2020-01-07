#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LEO1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LEO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LEP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LEPR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LEPROT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LEPROTL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LETM1 
