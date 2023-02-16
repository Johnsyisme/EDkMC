units metal
boundary p p p
atom_style atomic
neighbor 2.0 bin
neigh_modify every 1

variable j equal step_tmp

compute tmp all pe

read_data BN.data
pair_style deepmd ../1029c.pb
pair_coeff * *
thermo     1
thermo_style custom step pe temp c_tmp
minimize 1.0e-7 1.0e-10 1000 10000
write_data relaxed.data
run 0

variable E_new equal c_tmp
variable E python kmc
python kmc input 3 v_E_new 0.0 v_j return v_E format ffif file kMC_maBN.py
variable E_old equal ${E}

label      loop
variable   i loop 200000
variable   k equal v_i+step_tmp
clear

units metal
boundary p p p
atom_style atomic
neighbor 2.0 bin
neigh_modify every 1
compute tmp all pe
read_data BN.data
pair_style deepmd ../1029c.pb
pair_coeff * *
thermo     1
thermo_style custom step pe temp c_tmp
minimize 1.0e-7 1.0e-10 1000 10000
write_data relaxed.data
run 0
variable E python kmc
variable E_new equal c_tmp
python kmc input 3 v_E_new v_E_old v_k return v_E format ffif file kMC_maBN.py
variable E_old equal ${E}

next i 
jump SELF loop
label      break
quit
