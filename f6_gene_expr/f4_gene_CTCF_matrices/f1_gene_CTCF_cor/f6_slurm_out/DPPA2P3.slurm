#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DPPA2P3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DPPA2P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPPA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPPA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPPA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPRX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPRXP4 
