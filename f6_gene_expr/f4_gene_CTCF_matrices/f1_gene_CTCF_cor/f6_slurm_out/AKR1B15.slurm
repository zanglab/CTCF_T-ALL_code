#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AKR1B15.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AKR1B15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKR1C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKR1C2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKR1C3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKR1C4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKR1C6P 
