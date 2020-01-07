#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CTAGE9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CTAGE9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTBS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTCF 
