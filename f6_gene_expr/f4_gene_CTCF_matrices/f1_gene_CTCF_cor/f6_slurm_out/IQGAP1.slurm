#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IQGAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IQGAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQGAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQGAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQSEC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQSEC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQSEC3 
