#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LARP7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LARP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LARS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LARS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAS1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LASP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAT 
