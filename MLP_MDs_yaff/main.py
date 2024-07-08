from pathlib import Path
import parsl
from parsl.data_provider.files import File
import psiflow
import os
import molmod.units
import numpy as np
from ownscripts.bash_app_python import bash_app_python
from ownscripts.utils import run_MD

from utils_debug import plotvol, writeatoms, printenergy, writeinputtraj, create_ave_struc, get_bondangles_distribution

app_plotvol = bash_app_python(plotvol, precommand = "", executors=['default_threads'])
app_writeatoms = bash_app_python(writeatoms, precommand = "", executors=['default_threads'])
app_printenergy = bash_app_python(printenergy, precommand = "", executors=['default_threads'])
app_writeinputtraj = bash_app_python(writeinputtraj, precommand = "", executors=['default_threads'])
app_create_ave_struc = bash_app_python(create_ave_struc, precommand = "", executors=['default_threads'])
app_get_bondangles_distribution = bash_app_python(get_bondangles_distribution, precommand = "", executors=['default_threads'])
app_run_MD = bash_app_python(run_MD, precommand = "", executors=['ModelEvaluation'])

def app_run_MD_general(output_folder, simname, path_atoms, mult_steps = 1, seed = 0):
    #path calculator
    path_calc = File(str(data_folder / "OwnMLP.model"))
    #General settings MD
    total_steps = 200000 # for some simulation this number is doubled
    freq_step = 100
    temperature = 300
    pressure = 0.1*molmod.units.pascal*(10**6)
    
    path_output = File(str(output_folder / str("traj_"+simname+".xyz")))
    if os.path.exists(path_output.filepath) == False:
        MDapp = app_run_MD(
            execution_folder = output_folder,
            stderr           = str(output_folder / str("error_MD_"+simname+".txt")),
            stdout           = str(output_folder / str("output_MD_"+simname+".txt")),
            inputs           = [path_atoms, path_calc], 
            outputs          = [path_output], 
            steps            = int(mult_steps*total_steps), 
            step             = freq_step, 
            temperature      = temperature, 
            pressure         = pressure, 
            seed             = seed, 
            device           = "cuda")
        path_output = MDapp.outputs[0]
    return path_output

def app_create_ave_struc_general(output_folder, path_traj, simname):
    #General settings MD
    calib = 100 # number of printed steps (so multiply with freq_step for total steps) used as calibration and therefor not included in any post analysis
    path_ave_atoms = File(str(output_folder / str("ave_atoms_"+simname+".cif")))
    if os.path.exists(path_ave_atoms.filepath) == False:
        ave_at_app = app_create_ave_struc(
            execution_folder = output_folder,
            stderr           = str(output_folder / str("error_aveat_"+simname+".txt")),
            stdout           = str(output_folder / str("output_aveat_"+simname+".txt")),
            inputs           = [path_traj],
            outputs          = [path_ave_atoms],
            calib            = calib)
        return ave_at_app

def app_writeatoms_general(output_folder, path_traj, simname, index = -1):
    if index == -1:
        simname = "end_" + simname
    elif index == 0:
        simname = "begin_" + simname
    else:
        simname = "at"+str(index)+"_" + simname
    path_output = File(str(output_folder / str("atoms_"+simname+".xyz")))
    if os.path.exists(path_output.filepath) == False:
        wa_app = app_writeatoms(
            execution_folder = output_folder,
            stderr           = str(output_folder / str("error_write_"+simname+".txt")),
            stdout           = str(output_folder / str("output_write_"+simname+".txt")),
            inputs           = [path_traj], 
            outputs          = [path_output],
            index            = index)
        path_output = wa_app.outputs[0]
    return path_output


if __name__ == '__main__':
    psiflow.load()

    create_input_struc = True
    calib = 100 # number of printed steps (so multiply with freq_step for total steps) used as calibration and therefor not included in any post analysis

    #Create input and output folders
    Path_folder = Path.cwd() 
    data_folder = Path_folder / "data"
    main_output_folder = Path_folder / "output"
    main_output_folder.mkdir(exist_ok=True)
    print("Start constructing the simulations")

    labels_ba = []
    inputs_ba = []
    for ns in range(0,33,4):
        labels_ba.append(str(ns))
        for sd in range(10):
            simname = "OwnMLP_"+str(ns)+"_"+str(sd)
            path_atoms = File(str(data_folder / str("atoms_end_DMA_subst_"+str(ns)+"_"+str(sd)+".xyz")))
            outputtraj = app_run_MD_general(main_output_folder, simname, path_atoms, seed = sd)
            inputs_ba.append(outputtraj)
            app_create_ave_struc_general(main_output_folder, outputtraj, simname)
            if create_input_struc:
                app_writeatoms_general(main_output_folder, outputtraj, simname)

    #Get Pb-I-Pb bond angles distribution
    bondangle_plot_file = File(str(main_output_folder / str("Bondangles_comparison.pdf")))
    if os.path.exists(bondangle_plot_file.filepath) == False:
        ba_app_1s = app_get_bondangles_distribution(
            execution_folder = main_output_folder,
            stderr           = str(main_output_folder / str("error_ba.txt")),
            stdout           = str(main_output_folder / str("output_ba.txt")),
            inputs           = inputs_ba, 
            outputs          = [bondangle_plot_file], 
            labels           = labels_ba,  
            calib            = calib)

    print("get results of apps")
    parsl.wait_for_current_tasks()

