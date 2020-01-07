#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NINL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NINL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NIP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NIPA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NIPA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NIPAL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NIPAL2 
