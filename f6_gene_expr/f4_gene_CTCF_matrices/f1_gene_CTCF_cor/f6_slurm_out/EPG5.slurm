#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EPG5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EPG5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPGN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHA10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHA3 
