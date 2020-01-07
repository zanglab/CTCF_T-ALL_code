#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BTBD8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BTBD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTBD9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTF3L4 
