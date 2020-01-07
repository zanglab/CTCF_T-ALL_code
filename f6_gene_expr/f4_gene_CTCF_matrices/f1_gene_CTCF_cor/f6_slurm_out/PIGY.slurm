#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PIGY.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIH1D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIH1D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIH1D3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIK3AP1 
