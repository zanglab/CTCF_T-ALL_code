#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SFTA3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SFTA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFTPA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFTPA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFTPB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFTPC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFTPD 
