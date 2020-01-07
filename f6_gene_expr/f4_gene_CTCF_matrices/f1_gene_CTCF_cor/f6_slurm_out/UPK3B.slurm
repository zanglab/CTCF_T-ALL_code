#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UPK3B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UPK3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UPP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UPP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UPRT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UQCC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UQCC2 
