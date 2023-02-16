#!/bin/bash - 

module load lammps

mkdir structure
lmp -i in.lammps

a=` tail -1 ./energy_MC.log | awk '{print $2}'`
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
