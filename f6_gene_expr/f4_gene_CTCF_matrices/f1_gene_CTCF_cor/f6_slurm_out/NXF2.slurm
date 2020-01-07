#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NXF2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NXF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NXF2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NXF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NXF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NXF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NXN 
