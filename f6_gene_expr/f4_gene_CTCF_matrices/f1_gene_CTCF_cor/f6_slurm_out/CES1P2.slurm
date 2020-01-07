#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CES1P2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CES1P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CES2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CES3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CES4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CES5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CES5AP1 
