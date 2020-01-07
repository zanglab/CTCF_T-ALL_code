#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MDGA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MDGA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MDGA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MDH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MDH1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MDH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MDK 
