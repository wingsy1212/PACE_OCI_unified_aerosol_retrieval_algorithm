# PACE_OCI_unified_aerosol_retrieval_algorithm
Unified algorithm combining DT, DB, and OMI UV aerosol algorithm to be applied to PACE OCI sensor.

To run the program under run directory use run_PACE_Simulation.bash

Before run set up external data path:
under run/setup_env.bash  there needs to be directory for GMAO data

in run_PACE_Simulation.bash there needs to be change for related file names

line 77 path_obs  is the L1B input directory
line78 needs to change depending on the L1B file name

line94 call python depending on your system

line113/114 change ncgen directory
line 120 nccoy directory change


running the code is 
./run_PACE_Simulation.bash

Also under 
run/OCIUAAER_Config input output directory needs to be defined depending on the directory change of location of the dependency data, these data area shared using UMBC google drive within the release 
[https://drive.google.com/file/d/1T8na0ja0FDkkKAFsINZwTSs_pEivT2Br/view?usp=sharing](https://github.com/wingsy1212/PACE_OCI_unified_aerosol_retrieval_algorithm/releases)https://github.com/wingsy1212/PACE_OCI_unified_aerosol_retrieval_algorithm/releases![image](https://github.com/wingsy1212/PACE_OCI_unified_aerosol_retrieval_algorithm/assets/68506408/a07d73d0-66be-48a1-9954-6aa089549a99)


Under build folder needs to change Makefile based on the dependency.
