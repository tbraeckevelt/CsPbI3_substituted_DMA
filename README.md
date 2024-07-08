# CsPbI3_substituted_DMA
A repository for the paper about the effect of DMA substitutions on the structural properties of CsPbI3 

This repository consists of three folders:

- Create_init_struc: contains the files to generate initial structure to perform the active learning
*Create_subst.py: Python file to generate DMA (randomly rotated) substituted CsPbI3 with varying concentrations, uses GammaOpt.vasp and DMA_aligned.xyz as input
*UMLP_MACE.tar.gz: the files created with python script are a bit artificial and some relexation is necessary to generate realistic sample. 
                   To not perform to many acitve learning steps to perform this relexation, the foundation model of MACE (MACE-MP-0) is used instead (see MD simulation folder, 
                   the same MDs were performed with foundation model). The file is tarred te reduce memeory requirements.

- psiflow_files: contains the input files for the active learning via psiflow.
*hortense.yaml: contains the configuration settings for parsl, this selects the requested resources for the three different type of simulations (AI reference, MLP training, and MLP MD)
*run_sequential_learning.py: the main python script, which performs the psiflow active learning workflow.
*submitscript.sh: submits the python job to our local HPC, and loads the psiflow module (version 3.0.4)
*data: contains the input settings for CP2K, and the initial structures dataset (see also create_init_struc folder)
*output.txt: gives an overview of all MDs and generated stucture of the psiflow calculation
*output_final_iteration: contains the training and validation set, the log files, and the deployed and undeployed model (tarred) of the last iteration of the active learning loop

- MLP_MDs_yaff: contains the input files for the MDs performed with the final MLP constructed with psiflow. 
                A complete workflow is provided that performs the MD and calculates the Pb-I-Pb bond angle destributions.
*data: contains the initial structure for each MD (the same are used as for the psiflow active learning simulations). Moreover also the MACE model is given (tarred to reduce the size of the file).  
*main.py: python script that contains the full workflow (using parsl similar to psiflow active learning) to perform the MDs and to generate the bond angle destribution plot.
*utils_debug.py: contains some extra methods used in main.py.  
*hortense.yaml: contains the configuration settings for parsl, only the ModelEvaluation is used for the workflow.
*submitscript.sh: submits the python job to our local HPC, and loads and adapted psiflow module (version 3.0.0).  
*output_ba_000.txt and output_ba_001.txt: contains output of the generation of the bond angle plots, the average bond angle for different numbers of DMA substitutions is printed.
*Bondangles_comparison.pdf: plot showing the Pb-I-Pb bond angles destribution as a function of the number of DMA substituions.
