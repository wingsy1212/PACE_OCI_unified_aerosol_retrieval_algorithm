#!/bin/bash
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

function GetJulianDay # year, month, day
{
  year=$1
  month=$2
  day=$3

  jd=$((day - 32075 + 1461 * (year + 4800 - (14 - month) / 12) / 4 + 367 * (month - 2 + ((14 - month) / 12) * 12) / 12 - 3 * ((year + 4900 - (14 - month) / 12) / 100) / 4))

  echo $jd
}

# Deal with commandline options.
#------------------------------------------------------------------------------
# if [ $# -eq 0 ];  
# then
#   echo "Usage: `basename $0` (-ymdtgc) <ROOT_DIR>"
#   exit 1      
# fi  
root_dir='../DB_DT_Land/DB_Land'

while getopts ":y:m:d:t:g:c:o:i:" OPTION
do
  	case $OPTION in
  		y	)	opt_year=${OPTARG};;
#   		m	)	opt_month=${OPTARG};;
      m	)	opt_month=${OPTARG};;
  		d	)	opt_mday=${OPTARG};;  #day of year
  		t	) opt_time=${OPTARG};;
  		i	) opt_min=${OPTARG};;
  		g ) opt_granule=${OPTARG};;
  		c ) opt_collection=${OPTARG};;
  		o ) opt_dir=${OPTARG};;
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

if [ -n "$opt_month" ];
then
  month=$opt_month
  #`printf "%02d" $opt_month`
fi

if [ -n "$opt_mday" ];
then
  mday=$opt_mday
  #`printf "%02d" $opt_mday`
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


echo $year $month $mday $time $opt_min

rc=0

l1a_dir=../DB_DT_Land/DB_Land/input_data
if [ ! -d $l1a_dir ]; 
then 
	echo "L1A input directory does not exist.  Exiting."
	exit 1
fi
rm ${l1a_dir}/*
geo_datetime=${year}${month}${mday}
geo_doy0=$(date -d "$month / $mday / $year" '+%j')
geo_doy1=$(date -d "$month / $mday / $year" '+%j')

# echo l1a_dir: $l1a_dir

# for geo_file in `find -L ${l1a_dir} -maxdepth 1 -name "VNP03MOD.${file_regex}*.nc"`;
# do
#   echo "Current file: $geo_file"
#   file_base=`basename $geo_file .nc`
# 
#   geo_datetime=`echo $file_base | cut -d"." -f 2,3`
#   echo datetime:${geo_datetime}
#   geo_year=${geo_datetime:1:4}
#   geo_doy=${geo_datetime:5:3}
#   geo_doy=$((10#$geo_doy)) # force decimal (base 10) #remove zero  
#   geo_time=${geo_datetime:9:4} 
#   add_doy=$((geo_doy - 1))   
#   cald=`date -d "$add_doy days $geo_year-01-01" +"%Y%m%d"`
#   geo_year0=${cald:0:4}
#   geo_month0=${cald:4:2}
#   geo_day0=${cald:6:2}
# 
#   add_doy1=$((geo_doy))   
#   cald=`date -d "$add_doy1 days $geo_year-01-01" +"%Y%m%d"`
#   geo_year1=${cald:0:4}
#   geo_month1=${cald:4:2}
#   geo_day1=${cald:6:2}  
#   
#   geo_doy0=$(date -d "$geo_month0 / $geo_day0 / $geo_year0" '+%j')
#   geo_doy1=$(date -d "$geo_month1 / $geo_day1 / $geo_year1" '+%j')  
# 
#   echo year, month, mday, time, doy: $geo_year0 $geo_month0 $geo_day0 $geo_time $geo_doy0
# 
#   anc_year=$geo_year0
#   anc_month=$geo_month0
#   anc_day=$geo_day0
#   anc_hour=${geo_time:0:2}
#   anc_min=${geo_time:2:2}
# 
# # -- Get OMPS data
#   viirs_hour=$((10#$anc_hour))
#   viirs_min=$((10#$anc_min))
# 
#   viirs_juldate=$(GetJulianDay $anc_year 10#${anc_month} 10#${anc_day})
#   viirs_juldate=`echo "scale = 5; $viirs_juldate+$viirs_hour/24+$viirs_min/24/60" | bc`
#   #echo $viirs_juldate
# 
#   omps_files=`find -L ${l1a_dir} -maxdepth 1 -name "OMPS-NPP_NMEV-L1B-p000_v2.0_*.h5"`
#   omps_file_num=`find -L ${l1a_dir} -maxdepth 1 -name "OMPS-NPP_NMEV-L1B-p000_v2.0_*.h5"|wc -l`
#   echo "Number of OMPS files: $omps_file_num"
# 
#   n=0
#   for omps_file in $omps_files
#   do
#     #echo "Current file: $omps_file"
#     file_base=`basename $omps_file .h5`
#   
#     omps_datetime=`echo $file_base | cut -d"_" -f 4`
#     omps_year=${omps_datetime:0:4}
#     omps_month=${omps_datetime:5:2}
#     omps_day=${omps_datetime:7:2}   
#     omps_hour=${omps_datetime:10:2}
#     omps_min=${omps_datetime:12:2}
# 
#     omps_month=$((10#$omps_month))
#     omps_day=$((10#$omps_day))
#     omps_hour=$((10#$omps_hour))
#     omps_min=$((10#$omps_min))
#     #echo year, month, day, hour, min: $omps_year $omps_month $omps_day $omps_hour $omps_min
# 
#     omps_juldate_tmp=$(GetJulianDay $omps_year $omps_month $omps_day)
#     omps_juldate=`echo "scale = 5; $omps_juldate_tmp+$omps_hour/24+$omps_min/24/60" | bc`
# 
#     time_diff=`echo "scale = 5; $omps_juldate-$viirs_juldate" | bc`
#     time_diff=`echo "scale = 5; $time_diff*24*60" | bc`
#     time_diff=${time_diff%.*} # to integer
#     
#     if [ $time_diff -gt -80 -a $time_diff -lt 10 ];
#     then
#       omps_file_fin=$omps_file
#       time_diff_fin=$time_diff
#     fi 
#   done
  
  post_year=${year}
  #============for standard retrieval===========
  #===un-comment below for standard retrieval===  
  case "$time" in
  00|01|02)
    pre_model_run='00'
    post_model_run='03'
    pre_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${year}${month}${mday}_${pre_model_run}00.*.nc4
    post_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${year}${month}${mday}_${post_model_run}00.*.nc4
    ;;    
  03|04|05)
    pre_model_run='03'
    post_model_run='06'
    pre_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${year}${month}${mday}_${pre_model_run}00.*.nc4
    post_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${year}${month}${mday}_${post_model_run}00.*.nc4
    ;;     
  06|07|08)
    pre_model_run='06'
    post_model_run='09'
    pre_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${year}${month}${mday}_${pre_model_run}00.*.nc4
    post_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${year}${month}${mday}_${post_model_run}00.*.nc4
    ;;
  09|10|11)
    pre_model_run='09'
    post_model_run='12'
    pre_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${year}${month}${mday}_${pre_model_run}00.*.nc4
    post_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${year}${month}${mday}_${post_model_run}00.*.nc4
    ;;
  12|13|14)
    pre_model_run='12'
    post_model_run='15'
    pre_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${year}${month}${mday}_${pre_model_run}00.*.nc4
    post_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${year}${month}${mday}_${post_model_run}00.*.nc4
    ;;
  15|16|17)
    pre_model_run='15'
    post_model_run='18'
    pre_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${year}${month}${mday}_${pre_model_run}00.*.nc4
    post_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${year}${month}${mday}_${post_model_run}00.*.nc4
    ;; 
  18|19|20)
    pre_model_run='18'
    post_model_run='21'
    pre_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${year}${month}${mday}_${pre_model_run}00.*.nc4
    post_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${year}${month}${mday}_${post_model_run}00.*.nc4
    ;;         
  21|22|23)    
    pre_model_run='21'
    post_model_run='00'
    pre_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${year}${month}${mday}_${pre_model_run}00.*.nc4

#   -- increment the date by 1 day to pick up the correct 0z model data.
    echo ${year}${month}${mday} UTC + 1 day
    tmp_date=$(date -ud "${year}${month}${mday} UTC + 1 day" +%Y%m%d)
    post_year=${tmp_date:0:4}
    tmp_month=${tmp_date:4:2}
    tmp_day=${tmp_date:6:2}

    geo_doy1=$(date -d "$tmp_month / $tmp_day / $post_year" '+%j')  
	  post_anc_base=GEOS.*.*.inst3_2d_asm_Nx.GEOS*.${post_year}${tmp_month}${tmp_day}_${post_model_run}00.*.nc4
    ;;
    
  *)
    echo "Invalid time detected. Something went wrong! $gmtco_time"
    continue
    ;;
  esac

  pre_anc_dir="/tis/modaps/allData/1/DFPITI3NXASM/${year}/${geo_doy0}/"
  post_anc_dir="/tis/modaps/allData/1/DFPITI3NXASM/${post_year}/${geo_doy1}/"
#   pre_anc_dir="/data/wvkim/GEOS5/${year}/${geo_doy0}/"
#   post_anc_dir="/data/wvkim/GEOS5//${post_year}/${geo_doy1}/"


  pre_anc_file=`find -L ${pre_anc_dir} -maxdepth 2 -name "${pre_anc_base}"`
  post_anc_file=`find -L ${post_anc_dir} -maxdepth 2 -name "${post_anc_base}"`
  echo pre_anc_file: $pre_anc_file
  echo post_anc_file: $post_anc_file  

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
  pre_anc_base=`basename ${pre_anc_file}`
  post_anc_base=`basename ${post_anc_file}`
  if [ ! -e ${l1a_dir}/${pre_anc_base} ]
  then
    cp -s $pre_anc_file ${l1a_dir}/${pre_anc_base}
#     cp $pre_anc_file ${l1a_dir}/${pre_anc_base}
  fi
  if [ ! -e ${l1a_dir}/${post_anc_base} ]
  then
    cp -s $post_anc_file ${l1a_dir}/${post_anc_base}
#     cp $post_anc_file ${l1a_dir}/${post_anc_base}
  fi
  pre_anc_file=`find -L ${l1a_dir} -maxdepth 1 -name "${pre_anc_base}"`
  post_anc_file=`find -L ${l1a_dir} -maxdepth 1 -name "${post_anc_base}"`
   
  case "$month" in
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

  swir_vis_surf_coeffs_annual=${root_dir}/tbls/sfc_tbl/viirs_snpp_2.2_surface_database_annual_v2_20210515_0.1.hdf

#   creation_dt=`date +%Y%m%d%H%M%S`  #calendar date
  creation_dt=`date +%Y%j%H%M%S`    #day of year

  output_file="Pace_L2_Land_DB.${geo_datetime}.${time}${opt_min}"
#   output_file=$opt_dir/${output_file}'.'${collection}'.'${creation_dt}'.nc'
  output_file=$opt_dir/${output_file}'.nc'
  
  echo 'output:' $output_file
  
  tmp_dir=../DB_DT_Land/DB_Land/tmp
#   config_file=`mktemp -p ../DB_DT_Land/DB_Land/tmp -t VIIRS_DB.XXXXXXXX`|| exit 1
  rm ${tmp_dir}/PACE.config
  config_file=${tmp_dir}/PACE.config
  cp ../DB_DT_Land/DB_Land/supportMOD_PR04DB/config ${config_file}
  
#   echo config file: ${config_file}
#   cat ../DB_DT_Land/DB_Land/supportMOD_PR04DB/config >> $config_file
  
  echo "gdas1=$pre_anc_file"    >> $config_file
  echo "gdas2=$post_anc_file"    >> $config_file
  echo "surf_coeffs=$surf_coeffs" >> $config_file
  echo "modis_surf_db=$modis_surf_db" >> $config_file
  echo "viirs_surf_db=$viirs_surf_db" >> $config_file
  echo "swir_vis_surf_coeffs=$swir_vis_surf_coeffs" >> $config_file
  echo "swir_vis_surf_coeffs_annual=$swir_vis_surf_coeffs_annual" >> $config_file
  echo "output_file=$output_file" >> $config_file
  echo "year=$year" >> $config_file
  echo "month=$month" >> $config_file
  echo "day=$mday" >> $config_file
  echo "platform=VIIRS" >> $config_file

  echo "done : DB config" ${config_file}
exit $rc

# All done
