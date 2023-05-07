#!/bin/sh
#
# 	Author: 	Corey Bettenhausen
#				      Science System and Applications, Inc
#				      Goddard Space Flight Center
#				      301.614.5383
#				      corey.bettenhausen@ssaihq.com
# 
#   HISTORY:
#   Updated by W Kim  Feb/05/2019 Modified for NRT retrieval.
#                           For NRT retrieval, line 191 and 309 should be un-commented.
#                           Line 180 and 259 should be commented
#
#------------------------------------------------------------------------------

USAGE='Usage: run_deep_blue_nasa.sh [-ymdtgc]'


# Deal with commandline options.
#------------------------------------------------------------------------------
if [ $# -eq 0 ];  
then
  echo "Usage: `basename $0` (-ymdtgc) <ROOT_DIR>"
  exit 1      
fi  

while getopts ":y:d:t:g:c:" OPTION
do
  	case $OPTION in
  		y	)	opt_year=${OPTARG};;
#   		m	)	opt_month=${OPTARG};;
  		d	)	opt_mday=${OPTARG};;  #day of year
  		t	) opt_time=${OPTARG};;
  		g ) opt_granule=${OPTARG};;
  		c ) opt_collection=${OPTARG};;
  	  '?')  echo "$0: invalid option -$OPTARG"
  	        echo "Usage: $USAGE"
  	        exit 1
  	        ;;
  	esac
done
shift $((OPTIND - 1))

if [ -n "$opt_year" ];
then
  year=`printf "%04d" $opt_year`
fi

# if [ -n "$opt_month" ];
# then
#   month=`printf "%02d" $opt_month`
# fi

if [ -n "$opt_mday" ];
then
  mday=`printf "%03d" $opt_mday`
fi

if [ -n "$opt_time" ];
then
  time=$opt_time
fi

if [ -n "$opt_granule" ];
then
  granule=$opt_granule
fi

if [ -n "$opt_collection" ];
then
  collection=`printf "%03d" $opt_collection`
else
  collection="001"
fi 

root_dir=$1
if [ ! -d $root_dir ]; 
then 
	echo "Root directory does not exist.  Exiting."
	echo $rood_dir
	exit 1
fi

# Make deep_blue & l2_merge
#------------------------------------------------------------------------------
local_dir=$PWD
machine_name=`uname -n`
case $machine_name in
gs613-crl*)
  cd ../src
  make 
  ;;
tanami*)
  cd ../src
  make 
  ;;
sipssci*)
  cd ../src
  make 
  ;;
*)
  echo 'Skip compiling'
  ;;
esac
cd $local_dir

# Check that $L1A_DIR exists.
#------------------------------------------------------------------------------
l1a_dir=${root_dir}/input_data
if [ ! -d $l1a_dir ]; 
then 
	echo "L1A input directory does not exist.  Exiting."
	exit 1
fi
echo l1a_dir: $l1a_dir



# Build up our VIIRS file regex  
#------------------------------------------------------------------------------
file_regex='A'
if [ -n "$year" ]; 
then
  file_regex=${file_regex}${year}
else
  file_regex=${file_regex}'????'
fi

# if [ -n "$month" ];
# then
#   file_regex=${file_regex}${month}
# else
#   file_regex=${file_regex}'??' 
# fi

if [ -n "$mday" ];
then
  file_regex=${file_regex}${mday}
else
  file_regex=${file_regex}'???'
fi

file_regex=${file_regex}'.'

if [ -n "$time" ];
then
  file_regex=${file_regex}${time}
else
  file_regex=${file_regex}'????'
fi

# Override all options if -g is defined.
 if [ -n "$granule" ];
 then
   file_regex=$granule
 fi                                                                                                                   

echo ${file_regex}
n_files=`find -L ${l1a_dir} -maxdepth 1 -name "VJ103MOD.${file_regex}*.nc"|wc -l`
if [ $n_files -eq 0 ];
then
  echo No files found for ${l1a_dir}/VJ103MOD.${file_regex}*.nc.
  exit 1
fi

rc=0
for geo_file in `find -L ${l1a_dir} -maxdepth 1 -name "VJ103MOD.${file_regex}*.nc"`;
do
  echo "Current file: $geo_file"
  file_base=`basename $geo_file .nc`

  geo_datetime=`echo $file_base | cut -d"." -f 2,3`
  echo datetime:${geo_datetime}
  geo_year=${geo_datetime:1:4}
  geo_doy=${geo_datetime:5:3}
  geo_doy=$((10#$geo_doy)) # force decimal (base 10) #remove zero  
  geo_time=${geo_datetime:9:4} 
  add_doy=$((geo_doy - 1))
    
  cald=`date -d "$add_doy days $geo_year-01-01" +"%Y%m%d"`
  geo_month=${cald:4:2}
  geo_day=${cald:6:2}
  
  echo year, month, mday, time, doy: $geo_year $geo_month $geo_day $geo_time $geo_doy
# -- Get GDAS ancilliary files before and after current granule.
# --------------------------------------------------------------
  anc_year=$geo_year
  anc_month=$geo_month
  anc_day=$geo_day
  anc_hour=${geo_time:0:2}
  anc_min=${geo_time:2:2}

  #==============for ntr retrieval==============
  #====un-comment below for nrt retrieval=======
#     year1=${geo_year}
#     year2=${geo_year}
#     doy1=${geo_doy}
#     doy2=${geo_doy}
#     hr1=${anc_hour}
#     hr2=${anc_hour}
#     
#     one=1
#     if [ $anc_min -ge 30 ];
#     then    
#       let "hr2=10#${anc_hour} + ${one}"
#       if [ ${anc_hour} -eq 23 ];   
#       then 
#         hr2=00
#         let "doy2=10#${geo_doy} + ${one}"
#         
#         if [ ${geo_year} -eq 2020 -o ${geo_year} -eq 2024 -o ${geo_year} -eq 2028 ];
#         then 
#           if [ ${geo_doy} -eq 366 ];   
#           then
#             doy2=001
#             let "year2=10#${geo_year} + ${one}"
#           fi
#         else
#           if [ ${geo_doy} -eq 365 ];   
#           then
#             doy2=001
#             let "year2=10#${geo_year} + ${one}"
#           fi         
#         fi
#       fi  
#     else
#       let "hr1=10#${anc_hour} - ${one}"
#       if [ ${anc_hour} -eq 0 ];   
#       then 
#         hr1=23
#         let "doy1=10#${geo_doy} - ${one}"
#         if [ ${geo_doy} -eq 1 ];   
#         then 
#           doy1=365
#           let "year1=10#${geo_year} - ${one}"
#           if [ ${geo_year} -eq 2021 -o ${geo_year} -eq 2025 -o ${geo_year} -eq 2029 ];
#           then
#             doy1=366
#           fi
#         fi
#       fi  
#     fi
#     
#     ydh1=${year1}-${doy1}T${hr1}:00:00
#     nrtgeos=`sips_search -a ${ydh1} ${ydh1} DFPT1NXSLV|tail -n 1`
#     nrtgeos_trim=$(echo $nrtgeos | tr "/" "\n")
#     nrtgeos_arr=($nrtgeos_trim)
#     sz=${#nrtgeos_arr[@]}
#     nrtgeos1=${nrtgeos_arr[sz-1]}
#     
#     ydh2=${year2}-${doy2}T${hr2}:00:00
#     nrtgeos=`sips_search -a ${ydh2} ${ydh2} DFPT1NXSLV|tail -n 1`
#     nrtgeos_trim=$(echo $nrtgeos | tr "/" "\n")
#     nrtgeos_arr=($nrtgeos_trim)
#     sz=${#nrtgeos_arr[@]}
#     nrtgeos2=${nrtgeos_arr[sz-1]}
  #=============================================    
  #==============end ntr retrieval==============
  
# by default, assume the anc files are from the same date as the GMTCO files.
	machine_name=`uname -n`
	case $machine_name in
	gs613-crl*)
		pre_anc_dir="/data14/wvkim/GEOS5/${anc_year}/"
		post_anc_dir="/data14/wvkim/GEOS5/${anc_year}/"
		;;
	tanami*)
		pre_anc_dir="/data/wvkim/GEOS5/${anc_year}/"
		post_anc_dir="/data/wvkim/GEOS5/${anc_year}/"
		;;
	sipssci*)
    #============for standard retrieval===========
    #===un-comment below for standard retrieval===
    mkdir  ${root_dir}/input_data
#     sips_search -a ${geo_year}-${geo_doy} ${geo_year}-$((geo_doy+1)) DFPITI3NXASM | xargs -n 1 -I {} rsync {} .
    sips_search --download ${root_dir}/input_data -a ${geo_year}-${geo_doy} ${geo_year}-$((geo_doy+1)) DFPITI3NXASM

    #=============================================    
    #============end standard retrieval===========

    #==============for ntr retrieval==============
    #====un-comment below for nrt retrieval=======
#     cd input_data
#       sips_search -a ${ydh1} ${ydh1} DFPT1NXSLV|tail -n 1 | xargs -n 1 -I {} rsync {} .  
#       sips_search -a ${ydh2} ${ydh2} DFPT1NXSLV|tail -n 1 | xargs -n 1 -I {} rsync {} .    
#     cd ..
    #=============================================    
    #==============end ntr retrieval==============
    pre_anc_dir="../input_data/"
    post_anc_dir="../input_data/"   

		;;
	*)
		pre_anc_dir="../input_data/"
		post_anc_dir="../input_data/"
		;;
	esac
  #============for standard retrieval===========
  #===un-comment below for standard retrieval===  
  case "$anc_hour" in
  00|01|02)
    pre_model_run='00'
    post_model_run='03'
    pre_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${anc_year}${anc_month}${anc_day}_${pre_model_run}00.*.nc4
    post_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${anc_year}${anc_month}${anc_day}_${post_model_run}00.*.nc4
    ;;    
  03|04|05)
    pre_model_run='03'
    post_model_run='06'
    pre_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${anc_year}${anc_month}${anc_day}_${pre_model_run}00.*.nc4
    post_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${anc_year}${anc_month}${anc_day}_${post_model_run}00.*.nc4
    ;;     
  06|07|08)
    pre_model_run='06'
    post_model_run='09'
    pre_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${anc_year}${anc_month}${anc_day}_${pre_model_run}00.*.nc4
    post_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${anc_year}${anc_month}${anc_day}_${post_model_run}00.*.nc4
    ;;
  09|10|11)
    pre_model_run='09'
    post_model_run='12'
    pre_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${anc_year}${anc_month}${anc_day}_${pre_model_run}00.*.nc4
    post_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${anc_year}${anc_month}${anc_day}_${post_model_run}00.*.nc4
    ;;
  12|13|14)
    pre_model_run='12'
    post_model_run='15'
    pre_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*${anc_year}${anc_month}${anc_day}_${pre_model_run}00.*.nc4
    post_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${anc_year}${anc_month}${anc_day}_${post_model_run}00.*.nc4
    ;;
  15|16|17)
    pre_model_run='15'
    post_model_run='18'
    pre_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*${anc_year}${anc_month}${anc_day}_${pre_model_run}00.*.nc4
    post_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${anc_year}${anc_month}${anc_day}_${post_model_run}00.*.nc4
    ;; 
  18|19|20)
    pre_model_run='18'
    post_model_run='21'
    pre_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*${anc_year}${anc_month}${anc_day}_${pre_model_run}00.*.nc4
    post_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${anc_year}${anc_month}${anc_day}_${post_model_run}00.*.nc4
    ;;         
  21|22|23)    
    pre_model_run='21'
    post_model_run='00'
    pre_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${anc_year}${anc_month}${anc_day}_${pre_model_run}00.*.nc4

#   -- increment the date by 1 day to pick up the correct 0z model data.
    echo ${geo_year}${geo_month}${geo_day} UTC + 1 day
    tmp_date=$(date -ud "${geo_year}${geo_month}${geo_day} UTC + 1 day" +%Y%m%d)
    tmp_year=${tmp_date:0:4}
    tmp_month=${tmp_date:4:2}
    tmp_day=${tmp_date:6:2}

	  post_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${tmp_year}${tmp_month}${tmp_day}_${post_model_run}00.*.nc4
    ;;
    
  *)
    echo "Invalid time detected. Something went wrong! $gmtco_time"
    continue
    ;;
  esac

  pre_anc_file=`find -L ${pre_anc_dir} -maxdepth 2 -name "${pre_anc_base}"`
  post_anc_file=`find -L ${post_anc_dir} -maxdepth 2 -name "${post_anc_base}"`
  #=============================================    
  #============end standard retrieval===========

  #==============for ntr retrieval==============
  #====un-comment below for nrt retrieval=======
#     echo "NRT retrieval with geos5 FP file"
#     pre_anc_base=${nrtgeos1}
#     post_anc_base=${nrtgeos2}
#     pre_anc_file=`find -L ${pre_anc_dir} -maxdepth 2 -name "${pre_anc_base}"`
#     post_anc_file=`find -L ${post_anc_dir} -maxdepth 2 -name "${post_anc_base}"` 
  #=============================================    
  #==============end ntr retrieval==============    

  pre_anc_base=`basename ${pre_anc_file}`
  post_anc_base=`basename ${post_anc_file}`
  echo pre_anc_file: $pre_anc_file
  echo post_anc_file: $post_anc_file
  echo pre_anc_base: $pre_anc_base
  echo post_anc_base: $post_anc_base

  if [ -z $pre_anc_file ];
  then
    echo No ancillary data found prior to granule: : ${pre_anc_file}
    continue
  fi
  
  if [ -z $post_anc_file ];
  then
    echo No ancillary data found after granule: : ${post_anc_file}
    continue
  fi

# create link to ancillary files in input_data directory and reset
# $pre_anc_file and $post_anc_file.
  if [ ! -e ${l1a_dir}/${pre_anc_base} ]
  then
    cp -s $pre_anc_file ${l1a_dir}/${pre_anc_base}
  fi
  if [ ! -e ${l1a_dir}/${post_anc_base} ]
  then
    cp -s $post_anc_file ${l1a_dir}/${post_anc_base}
  fi
  pre_anc_file=`find -L ${l1a_dir} -maxdepth 1 -name "${pre_anc_base}"`
  post_anc_file=`find -L ${l1a_dir} -maxdepth 1 -name "${post_anc_base}"`
  echo pre: $pre_anc_file ${pre_anc_base}
  echo post: $post_anc_file ${post_anc_base}

  case "$geo_month" in
  12|01|02)
    surf_coeffs=${root_dir}/tbls/sfc_tbl/viirs_surfcoeffs_2ndpfit_winter_20170809.hdf
    modis_surf_db=${root_dir}/tbls/sfc_tbl/aqua_modis_surfdb_winter_20110913.hdf
    viirs_surf_db=${root_dir}/tbls/sfc_tbl/viirs_surfdb_winter_20170327.hdf
    swir_vis_surf_coeffs=${root_dir}/tbls/sfc_tbl/viirs_snpp_swir_vs_vis_surf_coeffs_winter_all_2nd_minimum_20170809.hdf  

    ;;
  03|04|05)
    surf_coeffs=${root_dir}/tbls/sfc_tbl/viirs_surfcoeffs_2ndpfit_spring_20170809.hdf
    modis_surf_db=${root_dir}/tbls/sfc_tbl/aqua_modis_surfdb_spring_20110913.hdf
    viirs_surf_db=${root_dir}/tbls/sfc_tbl/viirs_surfdb_spring_20170327.hdf
    swir_vis_surf_coeffs=${root_dir}/tbls/sfc_tbl/viirs_snpp_swir_vs_vis_surf_coeffs_spring_all_2nd_minimum_20170809.hdf

    ;;
  06|07|08)
    surf_coeffs=${root_dir}/tbls/sfc_tbl/viirs_surfcoeffs_2ndpfit_summer_20170809.hdf
    modis_surf_db=${root_dir}/tbls/sfc_tbl/aqua_modis_surfdb_summer_20110913.hdf
    viirs_surf_db=${root_dir}/tbls/sfc_tbl/viirs_surfdb_summer_20170327.hdf
    swir_vis_surf_coeffs=${root_dir}/tbls/sfc_tbl/viirs_snpp_swir_vs_vis_surf_coeffs_summer_all_2nd_minimum_20170809.hdf

    ;;
  09|10|11)
    surf_coeffs=${root_dir}/tbls/sfc_tbl/viirs_surfcoeffs_2ndpfit_fall_20170809.hdf
    modis_surf_db=${root_dir}/tbls/sfc_tbl/aqua_modis_surfdb_fall_20110913.hdf
    viirs_surf_db=${root_dir}/tbls/sfc_tbl/viirs_surfdb_fall_20170327.hdf
    swir_vis_surf_coeffs=${root_dir}/tbls/sfc_tbl/viirs_snpp_swir_vs_vis_surf_coeffs_fall_all_2nd_minimum_20170809.hdf
 
    ;;
  *)
    echo "Invalid month detected. Something went wrong!"
    continue
    ;;
  esac

  l1b_file=`find -L ${l1a_dir} -maxdepth 1 -name "VJ102MOD.${geo_datetime}.*.nc"`
  echo ${geo_datetime}
  if [ -z $l1b_file ];
  then
    echo No NASA L1B file found: VJ102MOD.${geo_datetime}_*.nc
    continue
  fi
  
#   creation_dt=`date +%Y%m%d%H%M%S`  #calendar date
  creation_dt=`date +%Y%j%H%M%S`    #day of year
  
  output_file="AERDB_L2_VIIRS_JPS1.${geo_datetime}"
  output_file=$root_dir/l2_data/${output_file}'.'${collection}'.'${creation_dt}'.nc'
  
  echo $output_file
  
  config_file=`mktemp -p $root_dir/tmp -t VIIRS_DB.XXXXXXXX`|| exit 1
  echo config file: ${config_file}
  cat $root_dir/supportMOD_PR04DB/config >> $config_file
  
  echo "geo_m=$geo_file"    >> $config_file
  echo "l1b_m=$l1b_file"    >> $config_file
  echo "gdas1=$pre_anc_file"    >> $config_file
  echo "gdas2=$post_anc_file"    >> $config_file
  echo "surf_coeffs=$surf_coeffs" >> $config_file
  echo "modis_surf_db=$modis_surf_db" >> $config_file
  echo "viirs_surf_db=$viirs_surf_db" >> $config_file
  echo "swir_vis_surf_coeffs=$swir_vis_surf_coeffs" >> $config_file
  echo "output_file=$output_file" >> $config_file
  echo "year=$geo_year" >> $config_file
  echo "month=$geo_month" >> $config_file
  echo "day=$geo_day" >> $config_file
  echo "platform=VIIRS" >> $config_file

  $root_dir/src/deep_blue $config_file
  echo config file: $config_file
  if [ -e $output_file ];
  then
    mv $config_file ${output_file}.config

#   -- compress output file with h5repack.
    compressed=`mktemp`
    h5repack -f GZIP=5  $output_file $compressed
    mv $compressed $output_file
  else
    echo "No output file. Retrieval bombed out."
    rc=1    
  fi

done

exit $rc

# All done
