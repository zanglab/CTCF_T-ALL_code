#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KDELR3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KDELR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDM1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDM1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDM2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDM2B 
