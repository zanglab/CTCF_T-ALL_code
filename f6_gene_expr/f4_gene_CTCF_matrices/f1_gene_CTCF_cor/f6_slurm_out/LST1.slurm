#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LST1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LST1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LTA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LTA4H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LTB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LTB4R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LTB4R2 
