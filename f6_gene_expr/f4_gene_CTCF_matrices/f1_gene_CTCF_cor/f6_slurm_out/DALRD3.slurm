#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DALRD3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DALRD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAND5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAOA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAP3 
