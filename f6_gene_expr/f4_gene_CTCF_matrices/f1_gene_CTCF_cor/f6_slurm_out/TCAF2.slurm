#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TCAF2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TCAF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCAF2P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCAIM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCAM1P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCEA1 
