#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ISM2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ISM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ISOC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ISOC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ISPD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IST1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ISX 
