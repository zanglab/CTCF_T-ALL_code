#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TAS2R41.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R41 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R50 
