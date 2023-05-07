#!/bin/sh
AERO_L1A_ROOT='../input_data/aeronet_colo/viirs'
L1A_ROOT='../input_data'
L2_ROOT='../l2_data'

if [ ! -e ${L2_ROOT}/aeronet_colo ]
then
	mkdir ${L2_ROOT}/aeronet_colo
fi

while getopts ":j:" OPTION
do
  	case $OPTION in
  		j ) opt_cpus=${OPTARG};;
  		'?')  echo "$0: invalid option -$OPTARG"
  	        echo "Usage: $USAGE"
  	        exit 1
  	        ;;
  	esac
done
shift $((OPTIND - 1))

# loop through all AERONET sites
for site in `ls -d ${AERO_L1A_ROOT}/*`
#for asite in mexico_city
#do
#for site in `ls -d ${AERO_L1A_ROOT}/${asite}`
do
  AERO_SITE=`basename ${site}`
  echo ${AERO_SITE}
  AERO_L2_DIR=${L2_ROOT}/aeronet_colo/${AERO_SITE}
  if [ -e ${AERO_L2_DIR} ]
  then
	  continue
  fi
  
  /usr/bin/find -L ${L1A_ROOT} -maxdepth 1 -name "*.nc" | xargs -n 1 rm -f
  if [ $? -ne 0 ] 
  then 
	  exit $?
  fi
  /usr/bin/find -L ${site} -maxdepth 1 -name "*.nc" | xargs -n 1 -I % cp -lv % ${L1A_ROOT}/
  if [ $? -ne 0 ] 
  then
    echo "Failed to copy L1A data. Stopping."
    exit $?
  fi

  if [ -z "$opt_cpus" ]
  then
    ./run_deep_blue_nasa.sh ..
  else
    ./run_deep_blue_parallel_nasa.sh -j $opt_cpus ..
  fi  
  
  AERO_L2_DIR=${L2_ROOT}/aeronet_colo/${AERO_SITE}
  mkdir ${AERO_L2_DIR}
  if [ $? -ne 0 ] 
  then 
	echo "Failed to make L2 output directory. Stopping."
	exit 1
  fi	
  /usr/bin/find -L ${L2_ROOT} -maxdepth 1 -name "*.nc" | xargs -n 1 -I % mv -v % ${AERO_L2_DIR}
  if [ $? -ne 0 ] 
  then
	    echo "Failed to move L2 data. Stopping."
	    exit $?
  fi
done
#done

