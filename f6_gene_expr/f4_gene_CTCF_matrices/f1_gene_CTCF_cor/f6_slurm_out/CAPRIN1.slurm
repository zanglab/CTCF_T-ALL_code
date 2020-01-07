#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CAPRIN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPRIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPRIN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPSL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPZA1 
