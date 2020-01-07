#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RIPPLY1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RIPPLY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIPPLY2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIPPLY3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RITA1 
