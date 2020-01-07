#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTGER2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGER2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGER3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGER4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGES 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGES2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGES3 
