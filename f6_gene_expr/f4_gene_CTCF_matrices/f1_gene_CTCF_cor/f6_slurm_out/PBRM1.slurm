#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PBRM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PBRM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PBX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PBX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PBX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PBX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PBXIP1 
