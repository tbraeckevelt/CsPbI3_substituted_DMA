import numpy as np
from ase.io import read, write
from ase.build.supercells import make_supercell
from scipy.spatial.transform import Rotation as R
from ase import Atom, Atoms



DMA_struc = read("DMA_aligned.xyz")

com_DMA = DMA_struc.get_center_of_mass().copy()
for at in DMA_struc:
    at.position -= com_DMA

gamma_struc = read("GammaOpt.vasp")
gamma_sup = make_supercell(gamma_struc, np.array([[2,2,0],[-2,2,0],[0,0,2]]))

for num_sub in range(4,29,8):
    for ns in range(10):
        supcel = gamma_sup.copy()

        rd_perm = np.random.permutation(64)[:num_sub]
        ind_Cs = -1
        sym_lst= []
        pos_lst= []
        for i,at in enumerate(supcel):
            if at.symbol == "Cs":
                ind_Cs +=1
                if ind_Cs in rd_perm:
                    copy_DMA_struc = DMA_struc.copy()
                    rot_mat =  R.random().as_matrix()
                    for atDMA in copy_DMA_struc:
                        atDMA.position = np.dot(rot_mat, atDMA.position)
                        atDMA.position += at.position
                        sym_lst.append(atDMA.symbol)
                        pos_lst.append(atDMA.position.copy())
                else:
                    sym_lst.append(at.symbol)
                    pos_lst.append(at.position.copy())
            else:
                sym_lst.append(at.symbol)
                pos_lst.append(at.position.copy())

        new_atoms = Atoms(symbols=sym_lst, positions=pos_lst, pbc=True, cell=supcel.cell)
        
        #get shortest distance between atoms
        for i, at in enumerate(new_atoms):
            for j, at2 in enumerate(new_atoms):
                if at.symbol in ["H"] and at2.symbol in ["I"]:
                    dist = new_atoms.get_distance(i, j, mic=True, vector=False)
                    if dist < 2.0:
                        #print("Short distance between atoms in ", num_sub, "substitutions", ns, "th structure")
                        print(num_sub, "substitutions", ns, "th structure", i, j, dist)

        name = "DMA_subst_"+str(num_sub)+"_"+str(ns) +".xyz"
        write(name, new_atoms)