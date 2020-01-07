#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PAXBP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PAXBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAXIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PBDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PBK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PBLD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PBOV1 
