#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACVR1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACVR1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACVR1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACVR2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACVR2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACVRL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACY1 
