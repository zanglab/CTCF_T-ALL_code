#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TNFAIP8L2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFAIP8L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFAIP8L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF10A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF10B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF10C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF10D 
