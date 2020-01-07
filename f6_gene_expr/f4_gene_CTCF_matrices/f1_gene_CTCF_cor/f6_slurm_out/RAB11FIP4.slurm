#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAB11FIP4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB11FIP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB11FIP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB15 
