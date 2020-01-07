#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BPGM.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BPGM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BPHL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BPI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BPIFA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BPIFA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BPIFA3 
