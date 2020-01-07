#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ENY2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ENY2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EOGT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EOMES 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EP300 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EP400 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPAS1 
