#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LIMS3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LIMS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIN28A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIN28B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIN37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIN52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIN54 
