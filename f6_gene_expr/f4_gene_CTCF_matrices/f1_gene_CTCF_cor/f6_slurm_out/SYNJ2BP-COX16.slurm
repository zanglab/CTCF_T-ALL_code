#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SYNJ2BP-COX16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNJ2BP-COX16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNPO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNPO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNPO2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNPR 
