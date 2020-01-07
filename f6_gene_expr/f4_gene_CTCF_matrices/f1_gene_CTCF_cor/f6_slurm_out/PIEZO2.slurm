#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PIEZO2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PIEZO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIFO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGBOS1 
