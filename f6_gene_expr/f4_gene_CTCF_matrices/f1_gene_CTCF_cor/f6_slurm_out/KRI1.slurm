#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRI1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRIT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT12 
