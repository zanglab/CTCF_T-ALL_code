#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SUPT20HL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SUPT20HL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUPT20HL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUPT3H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUPT4H1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUPT5H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUPT6H 
