#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UACA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UACA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UAP1L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBA3 
