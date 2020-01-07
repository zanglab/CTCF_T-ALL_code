#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NGF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NGF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NGFR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NGLY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NGRN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NHEJ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NHLH1 
