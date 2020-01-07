#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o METTL18.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL21A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL21C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL21EP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL23 
