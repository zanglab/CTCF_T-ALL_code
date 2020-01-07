#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRAF4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAFD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAIP 
