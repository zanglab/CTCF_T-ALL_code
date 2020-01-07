#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARMCX5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMCX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMCX6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARNT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARNT2 
