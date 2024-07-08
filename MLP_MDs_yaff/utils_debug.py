
def plotvol(inputs= [], outputs = [], el_lst = [], subst_lst = [], ns_lst = [], sd_lst = [], calib = 0):
    from matplotlib import pyplot as plt 
    import numpy as np
    from ase.io import read

    tel = 0
    vol_dct = {}
    for el in el_lst:
        vol_dct[el] = {}
        for subst in subst_lst:
            vol_dct[el][subst] = {}
            for ns in ns_lst:
                vol_dct[el][subst][ns] = {}
                vol_dct[el][subst][ns]["runs"] = []
                for sd in sd_lst:
                    traj = read(inputs[tel].filepath, index = str(calib) + ":")
                    tel += 1
                    ave_vol = 0.0
                    for atoms in traj:
                        ave_vol += atoms.get_volume()/len(traj)
                    vol_dct[el][subst][ns]["runs"].append(ave_vol)
                vol_dct[el][subst][ns]["ave"] = np.average(vol_dct[el][subst][ns]["runs"])
                vol_dct[el][subst][ns]["std"] = np.std(vol_dct[el][subst][ns]["runs"])

    # Plotting section
    for el in el_lst:
        for subst in subst_lst:
            x_values = list(vol_dct[el][subst].keys())
            y_values = [vol_dct[el][subst][ns]["ave"] for ns in ns_lst]
            y_errors = [vol_dct[el][subst][ns]["std"] for ns in ns_lst]
            plt.errorbar(x_values, y_values, yerr=y_errors, fmt='-o', label=f'{el}-{subst}')

    plt.xlabel('Number of substitutions')
    plt.ylabel('Average Volume')
    plt.title('Average Volume of CsPbI3 as a function of the number of substitutions')
    plt.legend()
    plt.savefig(outputs[0].filepath)


def writeatoms(inputs= [], outputs = [], index = "-1"):
    from ase.io import read, write

    atoms = read(inputs[0].filepath, index = index)
    write(outputs[0].filepath, atoms)


def printenergy(inputs= [], outputs = [], plane_lst = [], el_lst = [], subst_lst = [], sd_lst = [], calib = 0):
    import numpy as np
    from ase.io import read

    tel = 0
    en_dct = {}
    for plane in plane_lst:
        en_dct[plane] = {}
        for el in el_lst:
            en_dct[plane][el] = {}
            for subst in subst_lst:
                en_dct[plane][el][subst] = {}
                en_dct[plane][el][subst]["runs"] = []
                for sd in sd_lst:
                    traj = read(inputs[tel].filepath, index = str(calib) + ":")
                    tel += 1
                    ave_en = 0.0
                    for atoms in traj:
                        ave_en += atoms.get_potential_energy()/len(traj)
                    en_dct[plane][el][subst]["runs"].append(ave_en)
                en_dct[plane][el][subst]["ave"] = np.average(en_dct[plane][el][subst]["runs"])
                en_dct[plane][el][subst]["std"] = np.std(en_dct[plane][el][subst]["runs"])
    
    #printing section
    f_out = open(outputs[0].filepath, 'w')
    for plane in plane_lst:
        f_out.write("------" + plane + "\n")
        str_head = "el"
        for subst in subst_lst:
            str_head += "                                      " + str(subst)
        f_out.write(str_head + "\n")
        for el in el_lst:
            str_el = el 
            for subst in subst_lst:
                str_el += "    " + str(en_dct[plane][el][subst]["ave"]) + " (" + str(en_dct[plane][el][subst]["std"])+ ")"
            f_out.write(str_el + "\n")
    f_out.close()
    

def writeinputtraj(inputs= [], outputs = []):
    from ase.io import read, write
    traj = []
    for input in inputs:
        atoms= read(input.filepath)
        traj.append(atoms)
    write(outputs[0].filepath, traj)


def create_ave_struc(inputs = [], outputs = [], calib = 0):
    from ase.io import read, write
    from ase import Atoms

    traj = read(inputs[0].filepath, index = str(calib) + ":")

    ave_pos = traj[0].get_positions()/ len(traj)
    ave_cell = traj[0].get_cell()[:]/ len(traj)
    ave_energy = traj[0].info['energy']/ len(traj)
    if len(traj) > 1:
        for atoms in traj[1:]:
            ave_pos += atoms.get_positions()/ len(traj)
            ave_cell += atoms.get_cell()[:]/ len(traj)
            ave_energy += atoms.info['energy']/ len(traj)

    ave_atoms = Atoms(
                    numbers = atoms.get_atomic_numbers(),
                    positions = ave_pos,
                    cell = ave_cell,
                    pbc = True,
                    )
    ave_atoms.info['energy'] = ave_energy

    write(outputs[0].filepath, ave_atoms)

def get_bondangles_distribution(inputs = [], outputs = [], labels = [],  calib = 0):
    from ase.io import read
    import numpy as np
    import matplotlib.pyplot as plt

    if len(inputs) != len(labels):
        num_runs = int(len(inputs)/len(labels))
        assert len(inputs) == num_runs * len(labels), "the number of labels is not a devisor of number of inputs"
    else:
        num_runs = 1

    traj_lst = []
    for f in inputs:
        traj_lst.append(read(f.filepath, index = str(calib) + ":"))

    fig, ax = plt.subplots()

    for tel, traj in enumerate(traj_lst):
        if (tel%num_runs) == 0:
            angles_all = np.array([])
        All_I_angles = []
        for at in traj[0]:
            if at.symbol == "I":
                Pb_bonds = [at.index]
                for at2 in traj[0]:
                    if at2.symbol == "Pb":
                        dist = traj[0].get_distance(at.index, at2.index, mic=True, vector=False)
                        if dist < 4.5:
                            Pb_bonds.append(at2.index)
                assert len(Pb_bonds) <= 3, "Found too many Pb neighbours (i.e. "+str(len(Pb_bonds))+") for I atom "+str(at.index)+" for traj "+labels[tel]
                if len(Pb_bonds) == 3:
                    All_I_angles.append(Pb_bonds)   #Remove the I atoms with only 1 Pb bond

        angle_np = np.zeros((len(traj), len(All_I_angles)))
        for i, atoms in enumerate(traj):
            for j, Pb_bond in enumerate(All_I_angles):
                angle_np[i, j] = atoms.get_angle(Pb_bond[1], Pb_bond[0], Pb_bond[2], mic = True)

        for j, Pb_bond in enumerate(All_I_angles):
            angles_all = np.concatenate((angles_all, angle_np[:,j]))

        if (tel%num_runs) == num_runs -1:
            print(labels[int(np.floor(tel/num_runs))], np.average(angles_all))
            ax.violinplot(angles_all, np.arange(0.5*np.floor(tel/num_runs), 0.5*np.floor(tel/num_runs) + 0.1), showmeans=True, showextrema=False, showmedians=False)

    # set style for the axes
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(0, 0.5*len(labels), 0.5))
    ax.set_xticklabels(labels)
    ax.set_xlim(-0.25, 0.5*len(labels) - 0.25)

    ax.set_xlabel("Material")
    ax.set_ylabel('Angle [degree]')
    ax.set_title('Pb-I-Pb angles distribution')

    plt.subplots_adjust(bottom=0.15, wspace=0.05)
    plt.savefig(outputs[0].filepath,bbox_inches='tight')
    plt.close()
