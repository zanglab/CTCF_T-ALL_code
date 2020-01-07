#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EID2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EID2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EID2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EID3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF1AD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF1AX 
