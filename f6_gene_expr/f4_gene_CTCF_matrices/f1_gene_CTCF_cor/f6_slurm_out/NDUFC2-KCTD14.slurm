#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NDUFC2-KCTD14.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFC2-KCTD14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFS4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFS5 
