#!/bin/bash 
# 
#$ -cwd 
#$ -V 
#$ -j y 
#$ -S /bin/bash 
# 
 
R --vanilla <./code/mainScript_revised.R> mainScript_revised.Rout
