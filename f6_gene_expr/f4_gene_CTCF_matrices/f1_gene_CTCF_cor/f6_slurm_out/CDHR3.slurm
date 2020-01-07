#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDHR3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDHR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDHR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDHR5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDIPT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDIPT-AS1 
