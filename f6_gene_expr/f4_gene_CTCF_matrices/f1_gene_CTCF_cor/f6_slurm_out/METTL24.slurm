#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o METTL24.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL4 
