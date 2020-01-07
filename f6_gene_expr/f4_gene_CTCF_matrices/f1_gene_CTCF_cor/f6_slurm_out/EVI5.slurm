#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EVI5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EVI5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EVI5L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EVL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EVPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EVPLL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EVX1 
