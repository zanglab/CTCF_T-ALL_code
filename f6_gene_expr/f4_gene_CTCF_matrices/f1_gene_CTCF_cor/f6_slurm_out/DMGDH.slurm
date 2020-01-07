#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DMGDH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DMGDH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMKN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMPK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMRT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMRT2 
