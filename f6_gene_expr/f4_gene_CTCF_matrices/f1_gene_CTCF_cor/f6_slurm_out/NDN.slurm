#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NDN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NDN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDNF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDOR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDRG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDRG2 
