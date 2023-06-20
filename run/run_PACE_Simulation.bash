#!/bin/bash
########################################################################
# Puopose: 
# Generate aerosol products with the ABI netcdf files and
#      the retrieval model package (e.g., Dark Target )
# 
# Author:  
#       Zhaohui Zhang
#
# Usg: run_abi_v1.bash YYYYMNDD|YYYYJDAY(e.g., 20180831) HH MM 
#
########################################################################
#ulimit -s unlimited
 
ScriptDIR=`pwd`
source setup_env.bash 
########################################################################
# Set directories  
########################################################################
cd ..
PKG_HOME=$PWD/DT_Ocean 
input_dir=$PWD/L1_data
build_src=$PWD/build
run_dir=$PWD/run
OUT_DATA_HOME=$PWD/L2_data
write_dir=$PWD/Write_L2_data/pkg_root
dir_includ1=$PWD/DT_Common_src/inc
dir_includ2=$PWD/inc 
UTL_DIR=${PKG_HOME}/pkg_root 
RUNDIR=${run_dir} 
cd $run_dir   
  file_Flag="SYNTH" 
   echo $file_Flag  
#     for filename in  `cat File_list.txt` ; do
     for filename in  `cat Testing_File_list.txt` ; do 
     echo "the next file is $filename"
      PACE=`echo $filename | cut -d. -f1 | cut -c1-12` 
      yyyy=`echo $filename | cut -d. -f2 | cut -c1-4`
       mn=`echo $filename | cut -d. -f2 | cut -c5-6`
       dd=`echo $filename | cut -d. -f2 | cut -c7-8`
       hh=`echo $filename | cut -d. -f2 | cut -c10-11`
       mm=`echo $filename | cut -d. -f2 | cut -c12-13`    
#      echo $yyyy,$mn,$dd,$hh,$mm
# specify package components
if [ "${PKG_HOME}" == "" ]; then
  echo "Package home directory is not specified ... "
  exit 1
fi

CDATE=${UTL_DIR}/Routines/read_Anc/bin/cdate
WGRIB=${UTL_DIR}/Routines/read_Anc/bin/wgrib
if [ ! -f ${RUNDIR}/cdate ]; then
  cp $CDATE ${RUNDIR}/cdate
  chmod 755 ${RUNDIR}/cdate
fi 

# make loop lists
ddd_list=($mn$dd)
hrs_list=($hh)
min_list=($mm) 
set echo
rm OCIUAAER_Config
 
 for ddd in "${ddd_list[@]}" ; do
   for hh in "${hrs_list[@]}" ; do
     for mm in "${min_list[@]}"; do
       year=$yyyy
       month=${ddd:0:2} ; day=${ddd:2:2}
       jday=`date -d "$year$month$day" +%j`
       sdate=$year$jday
      
       echo "start" 
          echo $year,$month,$day,$hh$mm,$jday,$sdate,$ddd,$dd     
          rm $dir_includ1/read_Sat_MODIS.inc
          rm $dir_includ2/output_Variables.inc
#         path_obs=/data4/smattoo/Pace/Input/SyntheticData/L1B/$year/$jday 
          path_obs=$input_dir
          obsfile=${path_obs}/PACE_OCI_SIM.$year$month$day'T'$hh$mm*'.L1B.V9.1.nc' 
          echo input_l1file = $obsfile>input_obsfile 
        cat input_obsfile OCIUAAER_Config_main>OCIUAAER_Config  
#        echo $input_obsfile 
        cp $dir_includ1/read_Sat_MODIS_synth.inc   $dir_includ1/read_Sat_MODIS.inc
        cp $dir_includ2/output_Variables_synth.inc $dir_includ2/output_Variables.inc
      
      
 #  GMAO     FILE processing     
########################################################################
# processing ancillary data  Using /usr/local/bin/python3
########################################################################
       cd $RUNDIR ; echo ">>>> Processing ancillary data...." 
       # Get the right Names  for GDAS File based on time
       if [ -f ${TMP_DIR}/GDAS_Out ]; then  rm ${TMP_DIR}/GDAS_Out ; fi  
 #     Finding closest file name based on 3 hourly GMAO data   
    /usr/local/bin/python3 ${UTL_DIR}/Routines/read_Anc/ancillary_time_GMAO.py $Anc_Dir_GMAO $year$month$day $hh $mm 3  
########################################################################
# build DB config
########################################################################
      cd ${RUNDIR}
      echo 'start DB config'
      ./db_config.sh -y $year -m $month -d $day -t $hh -i $mm -o $OUT_DATA_HOME
           
########################################################################
# Build 
########################################################################
             cd ${build_src}
                 rm *.o read.exe  *.mod 
                  make -f Makefile_Pace_Aer_V2  
########################################################################
# Create Netcdf files and Run retrievals              
########################################################################
             
        cd  $RUNDIR
        /usr/bin/ncgen -o PACE_output.nc  ${write_dir}/CDL/PACE_L2.cdl
        /usr/bin/ncgen -o Interm_file.nc  ${write_dir}/CDL/Interm_L2.cdl 
        echo ">>>> Run Unified Retrieval algorithm ...." 
        ${build_src}/OCIUAAER.exe 
########################################################################
# NetCDF file is put into Outfile directory with right satellite name.
######################################################################## 
        /usr/bin/nccopy -d2 PACE_output.nc      ${OUT_DATA_HOME}/Pace_L2_Merged_$year$jday.$hh$mm.nc
        rm ${OUT_DATA_HOME}/Pace_L2_Land_DB.*
        rm *.nc input*    
     done #   mm
   done  #   hh
done   #  day
done   #  for all files

exit 0


