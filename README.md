# PACE_Unified_Algorithm
Unified algorithm combining DT, DB, and OMI UV aerosol algorithm to be applied to PACE OCI sensor.

To run the program under run directory use run_PACE_Simulation.bash

Before run set up external data path:
under run/setup_env.bash  there needs to be directory for GMAO data
Anc_package_path= to where you will put these binary files

in run_PACE_Simulation.bash there needs to be change for related file names

line 77 path_obs  is the L1B input directory
line78 needs to change depending on the L1B file name

line94 call python depending on your system

line113/114 change ncgen directory
line 120 nccoy directory change


running the code is 
./run_PACE_Simulation.bash

Under build folder needs to change Makefile based on the dependency.

