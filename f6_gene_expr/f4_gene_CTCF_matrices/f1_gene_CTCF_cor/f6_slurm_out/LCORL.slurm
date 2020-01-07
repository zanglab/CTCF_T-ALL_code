#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LCORL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LCORL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCTL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LDAH 
