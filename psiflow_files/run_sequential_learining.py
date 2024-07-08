import requests
import logging
from pathlib import Path
import numpy as np

from ase.io import read

import psiflow
from psiflow.learning import SequentialLearning, load_learning
from psiflow.models import MACEModel, MACEConfig
from psiflow.reference import CP2KReference
from psiflow.data import FlowAtoms, Dataset
from psiflow.walkers import DynamicWalker, PlumedBias
from psiflow.state import load_state        # necessary for restarting a run
from psiflow.metrics import Metrics


def get_reference():
    with open(Path.cwd() / 'data' / 'cp2k_input.txt', 'r') as f:
        cp2k_input = f.read()
    return CP2KReference(cp2k_input_str=cp2k_input)

def get_mace_model():
    config = MACEConfig()
    config.num_channels     = 16     #default
    config.max_L            = 1      #default
    config.r_max            = 7.0
    config.max_ell          = 3      #default
    config.correlation      = 4
    config.num_interactions = 2      #default
    config.energy_weight    = 10.0   #default
    config.forces_weight    = 1.0    #default
    config.lr               = 0.01   #default
    config.batch_size       = 2
    config.max_num_epochs   = 120
    config.patience         = 20
    config.ema              = True
    config.ema_decay        = 0.99   #default
    return MACEModel(config)

def main(path_output):
    assert not path_output.is_dir()
    reference = get_reference()     # CP2K; PBE-D3(BJ); TZVP
    model = get_mace_model()      # NequIP; default model
    trajectory = read(str(Path.cwd() / 'data' / 'input_traj.xyz'), index=':')

    can_traj = []
    for atoms in trajectory:
        flowatoms = FlowAtoms.from_atoms(atoms)
        flowatoms.canonical_orientation()  # transform into conventional lower-triangular box
        can_traj.append(flowatoms)

    model.add_atomic_energy("H", reference.compute_atomic_energy("H", box_size=6))
    model.add_atomic_energy("C", reference.compute_atomic_energy("C", box_size=6))
    model.add_atomic_energy("N", reference.compute_atomic_energy("N", box_size=6))
    model.add_atomic_energy("Cs", reference.compute_atomic_energy("Cs", box_size=6))
    model.add_atomic_energy("I", reference.compute_atomic_energy("I", box_size=6))
    model.add_atomic_energy("Pb", reference.compute_atomic_energy("Pb", box_size=6))

    # set learning parameters and do pretraining
    learning = SequentialLearning(
            path_output=path_output,
            niterations=10,
            train_valid_split=0.9,
            metrics= Metrics("Firsttry", "DMAdoping"),
            train_from_scratch=True,
            pretraining_nstates=90,
            pretraining_amplitude_pos=0.05,
            pretraining_amplitude_box=0.05,
            error_thresholds_for_reset=(10, 200),  # in meV/atom, meV/angstrom
            error_thresholds_for_discard=(30, 600), # in meV/atom, meV/angstrom
            temperature_ramp=(100, 600, 5),
            )
#look formation energies
    # construct walkers; biased MTD MD in this case
    walkers = DynamicWalker.multiply(
            90,
            data_start=Dataset(can_traj),
            timestep=0.5,
            steps=1000,
            step=100,
            start=0,
            temperature=600,
            pressure = 0.1, # NPT in MPa
            max_excess_temperature=1000,
            distance_threshold=0.7
            )
    data_train, data_valid = learning.run(
            model=model,
            reference=reference,
            walkers=walkers,
            )


def restart(path_output):
    reference = get_reference()
    learning  = load_learning(path_output)
    model, walkers, data_train, data_valid = load_state(path_output, '2')
    data_train, data_valid = learning.run(
            model=model,
            reference=reference,
            walkers=walkers,
            initial_data=data_train + data_valid,
            )


if __name__ == '__main__':
    psiflow.load()
    path_output = Path.cwd() / 'output' # stores learning results
    main(path_output)
    #restart(path_output)
