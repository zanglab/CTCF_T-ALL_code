#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDC7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC73 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDCA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDCA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDCA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDCA5 
