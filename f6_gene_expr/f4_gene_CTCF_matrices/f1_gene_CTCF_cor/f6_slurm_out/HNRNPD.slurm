#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HNRNPD.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPDL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPH3 
