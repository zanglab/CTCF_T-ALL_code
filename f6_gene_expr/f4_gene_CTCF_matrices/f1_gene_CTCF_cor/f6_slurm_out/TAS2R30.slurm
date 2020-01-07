#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TAS2R30.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R40 
