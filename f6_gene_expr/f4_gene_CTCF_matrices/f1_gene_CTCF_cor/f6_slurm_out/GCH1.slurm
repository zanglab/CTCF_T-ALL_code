#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GCH1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GCH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCHFR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCKR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCLC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCLM 
