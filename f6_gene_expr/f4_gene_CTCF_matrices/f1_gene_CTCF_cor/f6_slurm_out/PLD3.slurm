#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLD3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEK 
