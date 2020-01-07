#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLCB3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCE1 
