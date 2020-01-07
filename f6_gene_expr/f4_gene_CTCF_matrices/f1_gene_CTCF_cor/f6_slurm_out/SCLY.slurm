#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SCLY.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SCLY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCMH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCML1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCML2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCML4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCN10A 
