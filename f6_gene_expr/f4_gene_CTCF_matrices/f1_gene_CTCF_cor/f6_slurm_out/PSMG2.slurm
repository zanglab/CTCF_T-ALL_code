#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PSMG2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSORS1C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSORS1C2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSPC1 
