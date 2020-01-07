#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EIF4ENIF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF4ENIF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF4G1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF4G2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF4G3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF4H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF5 
