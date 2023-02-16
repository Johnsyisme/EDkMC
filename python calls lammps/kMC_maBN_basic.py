import numpy as np
import datetime
from pymatgen.core import Structure, Lattice
from pymatgen.io.vasp import Poscar
from scipy.spatial import cKDTree
import os

def readVasp(fn):
    return Poscar.from_file(fn).structure
def writeVasp(fn,structure):
    pos=Poscar(structure,comment='vasp')
    Poscar.write_file(pos,fn)
def writeFile(fn, s):
    with open(fn, 'w') as fout:
        fout.write(s)

def lastStructure():
    fn = "./accepted/step-%d.vasp"%i
    return readVasp(fn)
def randomStructure(supercell, lc):
    a, b = supercell
    area_prim = lc**2*3**0.5/2
    npairs = int((a*b)/area_prim)
    print("npairs: ", npairs)
    print("natoms: ", npairs*2)
    coords = np.random.rand(npairs*2, 3)
    coords[:, -1] = 0.5
    return Structure(lattice=[[a,0,0], [0,b,0], [0,0,20]], species=(["N"]*npairs+["B"]*npairs), \
        coords=coords, coords_are_cartesian=False)
def writeLammpsScript():
    out = [
        "units metal",
        "boundary p p p",
        "atom_style atomic",
        "read_data BN.data",
        "neighbor 2.0 bin",
        "neigh_modify every 1",
        "pair_style extep",
        "pair_coeff * * BN.extep N B",
        "thermo     1",
        "thermo_style custom step pe etotal",
        "minimize 1.0e-7 1.0e-10 1000 10000",
        "write_data relaxed.data"
    ]
    return writeFile("in.lammps", "\n".join(out))
def writeLammpsData(structure, fn):
    out = [
        "BN\n",
        "%d atoms\n"    % structure.num_sites,
        "%d atom types\n" % structure.ntypesp,
        "0.000000 %.6f  xlo xhi"%structure.lattice.a,
        "0.000000 %.6f  ylo yhi"%structure.lattice.b,
        "0.000000 %.6f  zlo zhi"%structure.lattice.c,
        "\nMasses\n",
        "1  14.0067",
        "2  10.8110", 
        "\nAtoms\n"
    ]
    for i,site in enumerate(structure.sites):
        out.append("%d %d %.6f %.6f %.6f"%((i+1, 1 if str(site.specie)=="N" else 2)  + tuple(list(site.coords))))
    writeFile(fn, "\n".join(out))
def data2Poscar(fn):
    with open(fn, 'r') as fin:
        data = fin.readlines()
        natoms = int(data[2].strip().split()[0])
        xmin, xmax = data[5].strip().split()[0],data[5].strip().split()[1]
        ymin, ymax = data[6].strip().split()[0],data[6].strip().split()[1]
        zmin, zmax = data[7].strip().split()[0],data[7].strip().split()[1]
    with open("relaxed.vasp",'w') as fout:
        fout.write("relaxed structure ma-BN\n")
        fout.write("1.0\n")
        fout.write("%f 0.000000 0.000000\n" % ( float(xmax)-float(xmin) ))
        fout.write("0.000000 %f 0.000000\n" % ( float(ymax)-float(ymin) ))
        fout.write("0.000000 0.000000 %f\n" % ( float(zmax)-float(zmin) ))
        fout.write("N  B\n")
        fout.write("%d  %d\n" % (int(0.5*natoms),int(0.5*natoms)))
        fout.write("Cartesian\n")
        for line in sorted(data[16:(16+natoms)],key=lambda line: line.split()[1]):     
            fout.write(" ".join(line.split()[2:5]) +"\n")
    return   
def parseLammps():
    data2Poscar("relaxed.data")
    with open("log.lammps", 'r') as fin:
        for line in fin:
            if "Energy initial, next-to-last, final =" in line: break
        line = fin.readline().strip().split()
        E_relaxed = float(line[2])
    return readVasp("relaxed.vasp"), E_relaxed
def runLammps():
    # lammps_exe = "srun lmp -in in.lammps"
    lammps_exe = r'D:\Program\"LAMMPS 64-bit 30Jul2021-MPI with Python"\bin\lmp.exe'
    #mpi_exe = r'D:\Program\MPICH2\bin\mpiexec.exe'
    os.system('echo 1')
    os.system("%s -in in.lammps 2>&1" %(lammps_exe) )
def findAtomPair(structure, cutoff=2.0):
    import itertools
    lattice = structure.lattice
    coords = structure.cart_coords
    pos = []
    for i, j in itertools.product([0,1],[0,1]):
        pos+=list(coords+[i*lattice.a, j*lattice.b, 0])
    kd = cKDTree(pos)
    pairs = np.array(list(kd.query_pairs(cutoff)))
    pairs = np.array(list(filter(lambda x: np.any(x<structure.num_sites), pairs)))
    pairs = pairs%structure.num_sites 
    return pairs[np.random.randint(len(pairs))]
def rotateSingleBond(structure):
    gindex1, gindex2 = findAtomPair(structure)
    pos1 = structure.cart_coords[gindex1]
    pos2 = structure.cart_coords[gindex2]
    center = (pos1+pos2)/2
    pos1 -= center
    pos2 -= center
    pos1 = [pos1[1], -pos1[0], 0] + center
    pos2 = [pos2[1], -pos2[0], 0] + center        
    structure.sites[gindex1].coords = pos1
    structure.sites[gindex2].coords = pos2
    return gindex1, gindex2
def exchangeBond(structure):
    gindex2=gindex1=0
    while structure.species[gindex1]==structure.species[gindex2]:
        gindex1, gindex2 = findAtomPair(structure)
    pos1, pos2 = structure.cart_coords[[gindex1, gindex2]]
    structure.sites[gindex2].coords=pos1
    structure.sites[gindex1].coords=pos2

def BNbondnum(structure):
    elementN1 = structure.cart_coords[list(structure.indices_from_symbol("N"))]
    elementB1 = structure.cart_coords[list(structure.indices_from_symbol("B"))]
    k_N1=cKDTree(elementN1)
    k_B1=cKDTree(elementB1)
    bondlength = 1.446 * 1.15
    index1 = k_N1.query_ball_tree(k_B1,bondlength)
    nBN1 = sum([len(i) for i in index1])
    return nBN1

def Run(nRot):
    os.system('mkdir unrelaxed accepted')
    structure = randomStructure([12.5621, 13.0550], 2.504)
    #structure = lastStructure()
    writeVasp("init.vasp", structure)
    writeLammpsScript()
    writeLammpsData(structure, "BN.data")
    runLammps()
    structure, E_new = parseLammps()
    writeVasp("relaxed-init.vasp", structure)
    kBT = 0.5
    E_old = E_new
    probability = np.random.rand(nRot,1)
    with open("energy_MC.log", "w+") as log:
        for nrot_ in range(nRot):
            oldStructure = structure.copy()
            nBN1=BNbondnum(oldStructure)
            if probability[nrot_] >= 0.5:
                rotateSingleBond(structure)
            else:
                exchangeBond(structure)
            writeVasp('./unrelaxed/step-%d.vasp'%(nrot_+i+1),structure)
            writeLammpsData(structure, "BN.data")
            runLammps()
            structure, E_new = parseLammps() 
            nBN2=BNbondnum(structure)
            E_BN= 2.15
            E_extra= (nBN2-nBN1)*E_BN         
            r = np.random.rand()
            p = np.exp(-(E_new-E_old-E_extra)/kBT) 
            s = ", ".join([
                    "Step: %6d"%(nrot_+i+1),
                    "Eold: %.4f"%E_old,
                    "Enew: %.4f"%E_new,
                    "E_diff: %.4f"%(E_new-E_old),
                    'nBN2: %4d'%nBN2 ,
                    'nBN2-nBN1: %4d' %(nBN2-nBN1),
                    "%s: %.3f"%("SW" if probability[nrot_] >= 0.5 else "EX",probability[nrot_]),
                    "%s:(%.3f, %.3f)" %("Denied" if r>p else "Accept",r, p)
                ])
            log.write(s+"\n")
            print(s)        
            log.flush() 
            if r>p:
                structure = oldStructure           
            else: 		
                writeVasp("./accepted/step-%d.vasp" % (nrot_+i+1),structure)     
                E_old = E_new   
    return                      
              
i = 0
starttime = datetime.datetime.now()
print("Programme start at: ",starttime)
if __name__ == "__main__":
	Run(200000)
endtime = datetime.datetime.now()
looptime = endtime - starttime
print ("CPU runtime: ", looptime)