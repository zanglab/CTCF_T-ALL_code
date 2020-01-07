#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IL17REL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IL17REL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL18BP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL18R1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL18RAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL19 
