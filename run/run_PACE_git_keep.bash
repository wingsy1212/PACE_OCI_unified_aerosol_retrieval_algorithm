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
input_dir=/tis/atmosphere/scratch/Dark_Target/smattoo/Projects_from_data4/PACE_github/Input 
input_dir=/tis/atmosphere/scratch/Dark_Target/smattoo/Projects_from_data4/PACE_unified_algorithm_from_paceteam/L1_data
OUT_DATA_HOME=$PWD/../Ouput
run_dir=$PWD/run  
Exec_dir=$PWD/bin 
echo $input_dir  
echo $OUT_DATA_HOME
echo $Anc_Dir_GMAO 
echo $Exec_dir
cd $run_dir
for filename in  `cat Testing_File_list.txt` ; do 
       echo "the next file is $filename"
      PACE=`echo $filename | cut -d. -f1 | cut -c1-12` 
      yyyy=`echo $filename | cut -d. -f2 | cut -c1-4`
       mn=`echo $filename | cut -d. -f2 | cut -c5-6`
       dd=`echo $filename | cut -d. -f2 | cut -c7-8`
       hh=`echo $filename | cut -d. -f2 | cut -c10-11`
       mm=`echo $filename | cut -d. -f2 | cut -c12-13`    
#      echo $yyyy,$mn,$dd,$hh,$mm
 # make loop lists
ddd_list=($mn$dd)
hrs_list=($hh)
min_list=($mm)  
 for ddd in "${ddd_list[@]}" ; do
   for hh in "${hrs_list[@]}" ; do
     for mm in "${min_list[@]}"; do
       year=$yyyy
       month=${ddd:0:2} ; day=${ddd:2:2}
       jday=`date -d "$year$month$day" +%j`
       sdate=$year$jday 
       echo "start" 
          echo $year,$month,$day,$hh$mm,$jday,$sdate,$ddd,$dd  
         path_obs=$input_dir
         obsfile=${path_obs}/PACE_OCI.$year$month$day'T'$hh$mm*'.L1B.nc'   
          echo   $obsfile                  
 #  GMAO     FILE processing     
########################################################################
# processing ancillary data  Using /usr/local/bin/python3
########################################################################
  echo ">>>> Processing ancillary data...."  
 # Finding closest file name based on 3 hourly GMAO data   
 python3 ancillary_time_2GMAO.py $Anc_Dir_GMAO $year$month$day $hh $mm 3 
 Anc1=`head -1 input_anc_file1`
 Anc2=`head -1 input_anc_file2` 
 
########################################################################
# Build 
########################################################################
 #            cd ${build_src}
 #                rm *.o read.exe  *.mod 
 #                 make -f Makefile_Pace_Aer_V2  
########################################################################
# Create Netcdf files and Run retrievals              
########################################################################
        /usr/bin/ncgen -o PACE_output.nc  PACE_L2.cdl
        /usr/bin/ncgen -o Interm_file.nc  Interm_L2.cdl 
        echo ">>>> Run Unified Retrieval algorithm ...." 
        ${Exec_dir}/oci_ua_aer ../uaa_lite.par $obsfile $Anc1 $Anc2 PACE_output.nc Interm_file.nc
########################################################################
# NetCDF file is put into Outfile directory with right satellite name.
######################################################################## 
    /usr/bin/nccopy -d2 PACE_output.nc   ${OUT_DATA_HOME}/Pace_L2_Merged_$year$jday.$hh$mm.nc
#        rm ${OUT_DATA_HOME}/Pace_L2_Land_DB.*
#         rm *.nc input*    
     done #   mm
   done  #   hh
done   #  day
done   #  for all files

exit 0


