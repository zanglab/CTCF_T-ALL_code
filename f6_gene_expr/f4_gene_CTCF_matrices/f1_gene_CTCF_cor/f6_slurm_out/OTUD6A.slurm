#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OTUD6A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OTUD6A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTUD6B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTUD7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTUD7B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTULIN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTX1 
