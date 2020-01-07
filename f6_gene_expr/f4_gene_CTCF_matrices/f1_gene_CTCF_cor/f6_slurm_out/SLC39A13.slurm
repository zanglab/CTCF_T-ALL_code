#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC39A13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC39A13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC39A14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC39A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC39A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC39A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC39A5 
