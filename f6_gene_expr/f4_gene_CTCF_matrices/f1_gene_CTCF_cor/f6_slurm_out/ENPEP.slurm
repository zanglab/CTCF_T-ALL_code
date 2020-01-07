#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ENPEP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ENPEP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENPP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENPP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENPP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENPP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENPP5 
