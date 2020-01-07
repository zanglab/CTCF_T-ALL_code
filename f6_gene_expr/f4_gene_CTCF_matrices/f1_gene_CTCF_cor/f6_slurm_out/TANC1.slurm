#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TANC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TANC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TANC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TANGO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TANGO6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TANK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAOK1 
