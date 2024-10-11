#!/bin/bash
#ulimit -s unlimited
  
ScriptDIR=`pwd`
source setup_env.bash 
########################################################################
# Set directories    
########################################################################
cd ..
#echo $PWD
dir="$PWD"
#echo $(dirname "$PWD") 
parentdir="$(dirname "$dir")"
#echo $parentdir
input_dir="$parentdir/Input"
echo $input_dir
input_dir=/run/cephfs/ACPS_Scratch/vkayetha/PACE-OCI/L1B
OUT_DATA_HOME=/tis/acps/scratch/vkayetha/PACE-OCI/GitVersion/Output
run_dir=$PWD/run  
Exec_dir=$PWD/bin 
echo $OUT_DATA_HOME
echo $Anc_Dir_GMAO 
echo $Exec_dir
cd $run_dir
#for filename in  `cat Testing_File_list.txt` ; do 
 for filename in  `cat File_list.txt` ; do 
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
         path_obs=$input_dir/$year/$month/$day/
         obsfile=$(find "${path_obs}" -name "PACE_OCI.${year}${month}${day}T${hh}${mm}*.L1B.V2.nc" )   
         echo   $obsfile 
 #  GMAO     FILE processing     
########################################################################
# processing ancillary data  Using /usr/local/bin/python3
########################################################################
  echo ">>>> Processing ancillary data...."  
 # Finding closest file name based on 3 hourly GMAO data   
/usr/local/bin/python3 ancillary_time_2GMAO.py $Anc_Dir_GMAO $year$month$day $hh $mm 3 
 Anc1=`head -1 input_anc_file1`
 Anc2=`head -1 input_anc_file2` 
 
########################################################################
# Build 
########################################################################
                  cd $Exec_dir
 #                rm *.o read.exe  *.mod 
                  make -f Makefile_Pace_Aer_V2 
		  cd $run_dir 
########################################################################
# Create Netcdf files and Run retrievals              
########################################################################
        rm   PACE_output.nc Interm_file.nc
        /usr/bin/ncgen -o PACE_output.nc  $PWD/../Write_L2_data/pkg_root/CDL/PACE_L2.cdl
        /usr/bin/ncgen -o Interm_file.nc  $PWD/../Write_L2_data/pkg_root/CDL/Interm_L2.cdl 
        echo ">>>> Run Unified Retrieval algorithm ...." 
        ${Exec_dir}/oci_ua_aer ../uaa_lite.par $obsfile $Anc1 $Anc2 PACE_output.nc Interm_file.nc
########################################################################
# NetCDF file is put into Outfile directory with right satellite name.
######################################################################## 
        /usr/bin/nccopy -d2 PACE_output.nc      ${OUT_DATA_HOME}/PACE_OCI.$year$month$day'T'$hh$mm.V2L1B.LatestChecks.nc
	
        rm *.nc    
     done #   mm
   done  #   hh
done   #  day
done   #  for all files

exit 0


