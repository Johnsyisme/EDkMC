#!/bin/bash - 
#===============================================================================
#
#          FILE: run.sh
# 
#         USAGE: ./run.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: YOUR NAME (), 
#  ORGANIZATION: 
#       CREATED: 09/22/2022 16:47
#      REVISION:  ---
#===============================================================================
module load lammps

mkdir structure
lmp -i in.lammps

a=` tail -1 ./energy_MC.log | awk '{print $2}'`
#a=` ls -l -r -t structure/step*| tail -1 | awk '{print $9}'`
c=${a%%\,*}

while [ $c -lt 200000 ]
do
	a=` tail -1 ./energy_MC.log | awk '{print $2}'`
	c=${a%%\,*}
	cp in.resume.bak in.resume
	cp bak.data BN.data
	sed -i "s/step_tmp/$c/g" in.resume
	lmp -i in.resume
done
