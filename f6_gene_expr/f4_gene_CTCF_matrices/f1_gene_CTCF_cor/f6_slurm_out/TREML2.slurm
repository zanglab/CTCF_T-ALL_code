#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TREML2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TREML2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TREML3P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TREML4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TREML5P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRERF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TREX1 
