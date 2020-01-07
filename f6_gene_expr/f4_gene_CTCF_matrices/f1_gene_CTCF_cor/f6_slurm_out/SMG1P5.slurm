#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMG1P5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMG1P5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMG1P7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMG5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMG6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMG7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMG8 
