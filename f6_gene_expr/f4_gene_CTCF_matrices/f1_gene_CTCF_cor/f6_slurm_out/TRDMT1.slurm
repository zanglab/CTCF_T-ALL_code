#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRDMT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRDMT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRDN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TREH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TREM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TREM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TREML1 
