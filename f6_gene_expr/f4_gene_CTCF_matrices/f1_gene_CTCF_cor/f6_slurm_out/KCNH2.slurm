#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCNH2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNH4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNH5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNH6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNH7 
