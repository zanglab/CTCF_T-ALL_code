#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHIT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHIT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHKA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHKB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHKB-CPT1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHM 
