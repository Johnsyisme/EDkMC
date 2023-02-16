import numpy as np
import datetime
from pymatgen.core import Structure, Lattice
from pymatgen.io.vasp import Poscar
from scipy.spatial import cKDTree
import os

def randomStructure(x,y):
    a, b = [x,y]
    area_prim = 2.504**2*3**0.5/2
    npairs = int((a*b)/area_prim)
    coords = np.random.rand(npairs*2, 3)
    coords[:, -1] = 0.5
    a=Structure(lattice=[[a,0,0],[0,b,0],[0,0,20]],species=(["B"]*npairs+["N"]*npairs) , \
        coords=coords, coords_are_cartesian=False)
    writeLammpsData(a,'BN.data')


def writeFile(fn, s):
    with open(fn, 'w') as fout:
        fout.write(s)
    return
def readFile(fn):
    with open(fn, 'r') as fin:
        return fin.read()
def readvasp(fn):
    return Poscar.from_file(fn).structure
def writevasp(structure,fn):
    pos=Poscar(structure)
    return Poscar.write_file(pos,fn)

def writeLammpsData(structure, fn):
    out = [
        "BN\n",
        "%d atoms"    % structure.num_sites,
        "%d atom types\n" % structure.ntypesp,
        "0.000000 %.6f  xlo xhi"%structure.lattice.a,
        "0.000000 %.6f  ylo yhi"%structure.lattice.b,
        "0.000000 %.6f  zlo zhi"%structure.lattice.c,
        "\nMasses\n",
        "1  10.8110",
        "2  14.0067",
        "\nAtoms\n"
    ]
    for i,site in enumerate(structure.sites):
        out.append("%d %d %.6f %.6f %.6f"%((i+1, 1 if str(site.specie)=="B" else 2) + tuple(list(site.coords))))
    writeFile(fn, "\n".join(out))

def parselammps(fn):
    with open(fn,'r') as fin2:
        lines=fin2.readlines()
        num=int(lines[2].strip().split()[0])
        [xmin,xmax]=lines[5].strip().split()[0:2]
        [ymin,ymax]=lines[6].strip().split()[0:2]
        [zmin,zmax]=lines[7].strip().split()[0:2]
        cord=np.array([i.strip().split() for i in lines[16:16+num]])
        specs=['B' if i=='1' else "N" for i in cord[:,1] ]
        coords=cord[:,2:5].astype(np.float64)
        lat=[[float(xmax)-float(xmin),0,0],[0,float(ymax)-float(ymin),0],[0,0,float(zmax)-float(zmin)]]
    a=Structure(lattice=lat,species=specs,coords=coords,coords_are_cartesian=True)
    a.sort()
    return a

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
    gindex1=0
    gindex2=0
    while structure.species[gindex1]==structure.species[gindex2]:
        gindex1,gindex2=findAtomPair(structure)
    pos1,pos2=structure.cart_coords[[gindex1,gindex2]]
    structure.sites[gindex2].coords=pos1
    structure.sites[gindex1].coords=pos2

def kmc(E_new,E_old,i):
    structure = parselammps('relaxed.data')
    kBT = 0.5
    log = open("energy_MC.log", "a+")
    r = np.random.rand(2,1)
    p = np.exp(-(E_new-E_old)/kBT)
    if r[0]>p:
        structure = parselammps('bak.data')
        s = ", ".join([
            "Step: %7d"%(i),
            "Eold: %12.5f"%E_old,
            "Enew: %12.5f"%E_new,
            "Denied (R, P): (%.3f, %.3f)" %(r[0], p)
        ])
    else:
        s = ", ".join([
            "Step: %7d"%(i),
            "Eold: %12.5f"%E_old,
            "Enew: %12.5f"%E_new,
            "Accepted (R, P): (%.3f, %.3f)" %(r[0], p)
        ])
        os.system('cp ./relaxed.data ./bak.data')
        E_old = E_new
        writevasp(structure,"./structure/step-%d.vasp"%(i))
    log.write(s+"\n")
    log.flush()
    log.close()
    if r[1]>=0.5:
        rotateSingleBond(structure)
    else:
        exchangeBond(structure)
    writeLammpsData(structure, "BN.data")
    return E_old
