#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NAA30.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NAA30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAA35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAA38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAA40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAA50 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAA60 
