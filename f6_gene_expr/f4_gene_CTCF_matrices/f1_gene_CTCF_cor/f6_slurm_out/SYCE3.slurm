#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SYCE3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SYCE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYCN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYCP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYCP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYCP2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYCP3 
