#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CWH43.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CWH43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CX3CL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CX3CR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXADR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXADRP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXADRP3 
