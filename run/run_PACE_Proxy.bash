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
build_src=$PWD/build
run_dir=$PWD/run
OUT_DATA_HOME=$PWD/L2_data
OUT_DATA_HOME1=$PWD/L2_data_forYingxi
write_dir=$PWD/Write_L2_data/pkg_root
dir_includ1=$PWD/DT_Common_src/inc
dir_includ2=$PWD/inc 
UTL_DIR=${PKG_HOME}/pkg_root  
RUNDIR=${run_dir} 
cd $run_dir 

if [ $# -ne  4 ]; then
 echo "usage: $0 YYYYMNDD|YYYYJDAY[e.g. 20181001|2018265] HH MM  PROXYorSYNTH"
 exit 1
fi
file_Flag=$(printf  "%05s"  $4)   
echo $file_Flag
 

if [ ${#1} -eq 8 ]; then
  PROC_TIME=`date +%Y%m%d  --date="$1"`
elif [ ${#1} -eq 7 ]; then
#  jdate="${1:0:4}0101 +$(( ${1:4:3} - 1 ))days"
jdate="${1:0:4}0101 +$((10#${1:4:3} - 1 ))days"
  PROC_TIME=$(date -d "$jdate" +%Y%m%d)
else
  echo "wrong date input format ...."
  exit 2
fi

yyyy=`echo ${PROC_TIME} | cut -c1-4`
mn=`echo ${PROC_TIME} | cut -c5-6`
dd=`echo ${PROC_TIME} | cut -c7-8`  
hh=$(printf "%02d" $((10#$2)))
mm=$(printf "%02d" $((10#$3)))
       echo "start one" 
      echo $yyyy,$mn,$dd,$hh,$mm
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
# Vinay Kayetha added option PROXY or SYNTHETIC.
        echo ${5}  
 for ddd in "${ddd_list[@]}" ; do
   for hh in "${hrs_list[@]}" ; do
     for mm in "${min_list[@]}"; do
       year=$yyyy
       month=${ddd:0:2} ; day=${ddd:2:2}
       jday=`date -d "$year$month$day" +%j`
       sdate=$year$jday
      
       echo "start"
    
      echo $year,$month,$day,$hh$mm,$jday,$sdate,$ddd,$dd   
      
      if [ $file_Flag == "PROXY" ]; then  
        rm $dir_includ1/read_Sat_MODIS.inc
        rm $dir_includ2/output_Variables.inc
 #       obsfile=TROP-in-Viirs_V4.1_A$sdate.$hh$mm'_TrORBIT-05582.h5'
 #       echo proxy_l1file = ../L1B_data/$obsfile>input_obsfile   
         path_obs=/data4/smattoo/Pace/Input/trop_in_viirs/$sdate 
 #        path_obs=/tis/acps/scratch/vkayetha/PACE_data_testbed/trop_in_viirs_fp/$sdate
         echo $path_obs
         obsfile=${path_obs}/TROP-in-Viirs_V4.1_A$sdate.$hh$mm*.h5
         echo $obsfile
         echo proxy_l1file = $obsfile>input_obsfile  
         echo $proxy_l1file 
        cat input_obsfile OCIUAAER_Config_main>OCIUAAER_Config  
         echo $input_obsfile
        cp $dir_includ1/read_Sat_MODIS_proxy.inc   $dir_includ1/read_Sat_MODIS.inc
        cp $dir_includ2/output_Variables_proxy.inc $dir_includ2/output_Variables.inc
        fi 
      if [ $file_Flag == "SYNTH" ]; then 
          rm $dir_includ1/read_Sat_MODIS.inc
          rm $dir_includ2/output_Variables.inc
         path_obs=/data4/smattoo/Pace/Input/SyntheticData/L1B/$year/$jday
         echo $path_obs
        obsfile=${path_obs}/PACE_OCI_SIM.$year$month$day'T'$hh$mm*'.L1B.V9.nc'
         echo $obsfile
        echo input_l1file = $obsfile>input_obsfile 
        cat input_obsfile OCIUAAER_Config_main>OCIUAAER_Config  
        echo $input_obsfile 
        cp $dir_includ1/read_Sat_MODIS_synth.inc   $dir_includ1/read_Sat_MODIS.inc
        cp $dir_includ2/output_Variables_synth.inc $dir_includ2/output_Variables.inc
      fi
      
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
# Input  input_anc_file  is created to be used by main programme.\
########################################################################
       echo ">>>> Processing GDS data...." $input_anc_file  
########################################################################
# build DB config
########################################################################
      cd ${RUNDIR}
      echo 'start DB config'
      ./db_config.sh -y $year -m $month -d $day -t $hh -i $mm -o $OUT_DATA_HOME
           
########################################################################
# Build & Run retrieval model 
########################################################################
             cd ${build_src}
                rm *.o read.exe  *.mod 
                 make -f Makefile_Pace_Aer_V2
             #exit 0
        cd  $RUNDIR
       /usr/bin/ncgen -o vnpaerdt_output_Land.nc ${UTL_DIR}/CDL/AERDT_L2_PACE_Land.cdl
       /usr/bin/ncgen -o vnpaerdt_output_Ocean.nc   ${UTL_DIR}/CDL/AERDT_L2_PACE_Ocean.cdl 
       /usr/bin/ncgen -o PACE_output.nc  ${write_dir}/CDL/PACE_L2.cdl
       /usr/bin/ncgen -o Interm_file.nc  ${write_dir}/CDL/Interm_L2.cdl
      
        echo ">>>> Run Retrieval model ...." 
        ${build_src}/OCIUAAER.exe 
########################################################################
# NetCDF file is put into Outfile directory with right satellite name.
######################################################################## 
 #    /usr/bin/nccopy -d2 vnpaerdt_output_Land.nc   ${OUT_DATA_HOME}/Pace_L2_Land_$year$jday.$hh$mm.nc
 #     /usr/bin/nccopy -d2 vnpaerdt_output_Ocean.nc  ${OUT_DATA_HOME}/Pace_L2_Ocean_$year$jday.$hh$mm.nc 
     /usr/bin/nccopy -d2 PACE_output.nc            ${OUT_DATA_HOME}/Pace_L2_Merged_$year$jday.$hh$mm.nc
#    /usr/bin/nccopy -d2   Interm_file.nc           $OUT_DATA_HOME1/Interm_L2_$year$jday.$hh$mm.nc 
#       rm *.nc input_* 
#        rm -r Temp_directory 
       
     done #   mm
   done  #   hh
done   #  day

exit 0


