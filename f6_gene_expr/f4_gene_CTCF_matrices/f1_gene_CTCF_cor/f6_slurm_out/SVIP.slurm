#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SVIP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SVIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SVOP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SVOPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SWAP70 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SWI5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SWSAP1 
