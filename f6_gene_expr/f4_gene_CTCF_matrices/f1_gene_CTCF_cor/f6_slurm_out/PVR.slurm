#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PVR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PVR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PVRIG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PVRIG2P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PWP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PWP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PWWP2A 
