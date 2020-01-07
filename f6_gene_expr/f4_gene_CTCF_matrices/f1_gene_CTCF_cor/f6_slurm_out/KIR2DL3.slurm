#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIR2DL3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIR2DL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIR2DL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIR2DS4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIR3DL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIR3DL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIR3DL3 
