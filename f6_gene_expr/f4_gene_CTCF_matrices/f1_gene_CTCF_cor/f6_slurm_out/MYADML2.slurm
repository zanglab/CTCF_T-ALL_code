#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYADML2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYADML2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYBBP1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYBL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYBL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYBPC1 
