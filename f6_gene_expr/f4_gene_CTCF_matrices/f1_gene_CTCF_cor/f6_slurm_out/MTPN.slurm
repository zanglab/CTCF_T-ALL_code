#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MTPN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MTPN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTRF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTRF1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTRNR2L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTRNR2L10 
