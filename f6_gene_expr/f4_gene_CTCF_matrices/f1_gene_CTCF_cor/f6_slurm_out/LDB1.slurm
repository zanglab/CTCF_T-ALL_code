#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LDB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LDB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LDB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LDB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LDHA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LDHAL6A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LDHAL6B 
