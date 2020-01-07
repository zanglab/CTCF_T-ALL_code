#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EIF4E1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF4E1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF4E2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF4E3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF4EBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF4EBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF4EBP3 
