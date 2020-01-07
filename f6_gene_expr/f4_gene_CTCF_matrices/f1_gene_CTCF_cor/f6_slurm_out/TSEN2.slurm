#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TSEN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TSEN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSEN34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSEN54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSFM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSG101 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSGA10 
