#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADPRHL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADPRHL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADPRM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADRA1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADRA1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADRA1D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADRA2A 
