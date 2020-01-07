#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DMRT3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DMRT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMRTA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMRTA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMRTB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMRTC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMRTC1B 
