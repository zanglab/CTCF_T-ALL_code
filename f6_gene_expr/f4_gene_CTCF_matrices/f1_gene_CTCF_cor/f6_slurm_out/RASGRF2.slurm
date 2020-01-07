#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RASGRF2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RASGRF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASGRP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASGRP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASGRP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASGRP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASIP1 
