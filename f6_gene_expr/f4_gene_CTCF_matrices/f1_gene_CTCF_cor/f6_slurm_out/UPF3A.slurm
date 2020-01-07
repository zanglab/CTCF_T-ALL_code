#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UPF3A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UPF3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UPF3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UPK1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UPK1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UPK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UPK3A 
