#!/bin/sh
#
# 	Author: 	Corey Bettenhausen
#				      Science System and Applications, Inc
#				      Goddard Space Flight Center
#				      corey.bettenhausen@ssaihq.com
#	
#    01 AUG 2017       JLee      Implemented machine-dependent path for parallel
#------------------------------------------------------------------------------

USAGE='Usage: run_deep_blue.sh [-cspeoq]'


# Deal with commandline options.
#------------------------------------------------------------------------------
if [ $# -eq 0 ];  
then
  echo "Usage: `basename $0` (-ymdtg) <ROOT_DIR>"
  exit 1      
fi  

while getopts ":y:m:d:t:g:j:" OPTION
do
  	case $OPTION in
  		y	)	opt_year=${OPTARG};;
#   		m	)	opt_month=${OPTARG};;
  		d	)	opt_mday=${OPTARG};;  #day of year
  		t	) opt_time=${OPTARG};;
  		g ) opt_granule=${OPTARG};;
  		j ) opt_cpus=${OPTARG};;
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

ncpus=1
if [ -n "$opt_cpus" ];
then
  ncpus=$opt_cpus
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
gs613-lakechad*)
  cd ../src
  make deep_blue_lake
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

n_files=`find -L ${l1a_dir} -maxdepth 1 -name "VNP03MOD.${file_regex}*.nc"|wc -l`
if [ $n_files -eq 0 ];
then
  echo No files found for ${l1a_dir}/VNP03MOD.${file_regex}*.
  exit 1
fi

dt=`find -L ${l1a_dir} -maxdepth 1 -name "VNP03MOD.${file_regex}*.nc" -exec basename {} \;|cut -d"." -f 2,3`
da=`find -L ${l1a_dir} -maxdepth 1 -name "VNP03MOD.${file_regex}*.nc" -exec basename {} \;|cut -d"." -f 2|cut -d"." -f 1`
uniq=($(printf "%s\n" "${da[@]}" | sort -u)); echo "${uniq[@]}"
nday=${#uniq[@]}
uniqdays=${uniq[@]}

case $machine_name in
  sipssci*)
  echo 'start GEOS5 sips_search'
  /home/wvkim/local/parallel-20170822/bin/parallel -j $ncpus ./geos5_sips_search.sh -d {} $root_dir ::: ${uniqdays[*]}  
  echo 'done GEOS5 sips_search'
  ;;
esac

machine_name=`uname -n`
case $machine_name in
gs613-crl*)
  /home/asayer/local/parallel-20170322/bin/parallel --delay 1.0 -j $ncpus ./run_deep_blue_nasa.sh -g {} $root_dir ::: ${dt[*]}
  ;;
tanami*)
  /home/cbettenh/local/parallel-20160622/bin/parallel --delay 1.0 -j $ncpus ./run_deep_blue_nasa.sh -g {} $root_dir ::: ${dt[*]}
  ;;
windhoek*)
  /home/wvkim/local/parallel/bin/parallel --delay 1.0 -j $ncpus ./run_deep_blue_nasa.sh -g {} $root_dir ::: ${dt[*]}
  ;;  
sipssci*)
  /home/wvkim/local/parallel-20170822/bin/parallel --delay 1.0 -j $ncpus ./run_deep_blue_nasa.sh -g {} $root_dir ::: ${dt[*]}  
  ;;
gs613-lakechad*)
  /home/wvkim/local_lib/parallel-20201122/bin/parallel --delay 1.0 -j $ncpus ./run_deep_blue_nasa.sh -g {} $root_dir ::: ${dt[*]}  
  ;;
*)
  echo "Invalid machine detected. Something went wrong!"
  continue
  ;;
esac


exit 1

# All done
