#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STRN3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STRN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STRN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STT3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STT3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STUB1 
