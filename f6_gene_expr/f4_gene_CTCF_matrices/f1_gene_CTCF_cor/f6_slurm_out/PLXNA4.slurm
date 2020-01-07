#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLXNA4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLXNA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLXNB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLXNB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLXNB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLXNC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLXND1 
