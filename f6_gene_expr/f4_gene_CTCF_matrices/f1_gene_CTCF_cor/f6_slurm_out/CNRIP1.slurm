#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CNRIP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CNRIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNST 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNTD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNTD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNTF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNTFR 
