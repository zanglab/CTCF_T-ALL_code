#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POLR3K.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR3K 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLRMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POM121 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POM121C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POM121L12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POM121L1P 
