#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CWC15.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CWC15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CWC22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CWC25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CWC27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CWF19L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CWF19L2 
