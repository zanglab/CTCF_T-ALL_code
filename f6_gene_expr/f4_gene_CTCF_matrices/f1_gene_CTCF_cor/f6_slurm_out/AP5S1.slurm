#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AP5S1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AP5S1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP5Z1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APAF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APBA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APBA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APBA3 
