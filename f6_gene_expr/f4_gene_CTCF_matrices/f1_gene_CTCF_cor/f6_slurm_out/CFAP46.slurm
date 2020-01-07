#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CFAP46.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP53 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP57 
