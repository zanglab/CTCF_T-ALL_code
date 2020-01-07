#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ISY1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ISY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ISY1-RAB43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ISYNA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITCH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITFG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITFG2 
