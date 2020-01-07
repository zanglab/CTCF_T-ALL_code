#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRAF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAF3IP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAF3IP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAF3IP3 
