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
 input_dir=/tis/atmosphere/scratch/Dark_Target/smattoo/Projects_from_data4/PACE_Git_updated_Sep2024/Input_L1
# input_dir=/tis/acps/scratch/vkayetha/PACE-OCI/L1B
OUT_dir=/tis/atmosphere/scratch/Dark_Target/smattoo/Projects_from_data4/PACE_Git_updated_Sep2024/V3_Output
#OUT_dir=/tis/atmosphere/scratch/Dark_Target/smattoo/Projects_from_data4/PACE_Git_updated_Sep2024/Output
echo $input_dir 

if [ ! -d ${OUT_dir} ];then
   mkdir -p ${OUT_dir}
fi


 
  
#OUT_DATA_HOME=$PWD/../Ouput
run_dir=$PWD/run  
Exec_dir=$PWD/bin 
echo $OUT_DATA_HOME
echo $Anc_Dir_GMAO 
echo $Exec_dir
cd $run_dir
   for filename in  `cat file_list_march23` ; do 
#     for filename in  `cat File_list.txt` ; do 
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
#           obsfile=$(find "${path_obs}" -name "PACE_OCI.${year}${month}${day}T${hh}${mm}*.L1B.V2.nc" )  
            obsfile=$(find "${path_obs}" -name "PACE_OCI.${year}${month}${day}T${hh}${mm}*.L1B_OCIT05.V3.nc" )   
           echo   $obsfile 
  OUT_DATA_HOME=$OUT_dir/$year/$month/$day/
 
   if [ ! -d ${OUT_DATA_HOME} ];then
   mkdir -p ${OUT_DATA_HOME}
fi
 echo $OUT_DATA_HOME
      
 #  GMAO     FILE processing     
########################################################################
# processing ancillary data  Using /usr/local/bin/python3
########################################################################
  echo ">>>> Processing ancillary data...."  
 # Finding closest file name based on 3 hourly GMAO data  
 rm input_anc_file1 input_anc_file2 
/usr/local/bin/python3 ancillary_time_GOESIT.py $Anc_Dir_GMAO $year$month$day $hh $mm 1 
Anc1=`head -1 input_anc_file1`
Anc2=`head -1 input_anc_file2` 
 echo $Anc1, $Anc2 
########################################################################
# Build 
########################################################################
 #            cd ${build_src}
 #                rm *.o read.exe  *.mod 
 #                 make -f Makefile_Pace_Aer_V2  
########################################################################
# Create Netcdf files and Run retrievals              
########################################################################
     rm   PACE_output.nc Interm_file.nc PACE_output_1KM.nc  
        /usr/bin/ncgen -o PACE_output.nc  PACE_L2.cdl
        /usr/bin/ncgen -o Interm_file.nc  Interm_L2.cdl 
        /usr/bin/ncgen -o PACE_output_1KM.nc PACE_L2_1KM.cdl 
        echo ">>>> Run Unified Retrieval algorithm ...." 
    echo $obsfile $Anc1,$Anc2, PACE_output.nc Interm_file.nc PACE_output_1KM.nc  
${Exec_dir}/oci_ua_aer ../uaa_lite.par $obsfile $Anc1 $Anc2 PACE_output.nc Interm_file.nc PACE_output_1KM.nc
########################################################################
# NetCDF file is put into Outfile directory with right satellite name.
######################################################################## 
        /usr/bin/nccopy -d2 PACE_output.nc      ${OUT_DATA_HOME}/Pace_L2_Merged_$year$jday.$hh$mm.nc
        /usr/bin/nccopy -d2 PACE_output_1KM.nc  ${OUT_DATA_HOME}/Pace_L2_Merged_1KM_$year$jday.$hh$mm.nc
 #       /usr/bin/nccopy -d2 Interm_file.nc      ${OUT_DATA_HOME}/Interm_$year$jday.$hh$mm.nc

#        rm *.nc    
     done #   mm
   done  #   hh
done   #  day
done   #  for all files

exit 0


