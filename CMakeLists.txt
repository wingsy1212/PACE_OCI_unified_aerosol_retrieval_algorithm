cmake_minimum_required(VERSION 3.12)

# set the Mac OSX target.  Normally not needed, but this
# is the only way I got the upgraded Monterey Mac to work
set(CMAKE_OSX_SYSROOT macosx)
set(CMAKE_OSX_DEPLOYMENT_TARGET "11")

project (PACE_UNIFIED_AEROSOL)

# where to look for the external libraries
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/opt/)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

message("CMake module path:" ${CMAKE_MODULE_PATH})

# set this to Release or Debug
# note: rerun cmake if you change the environment var OCSSW_DEBUG
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)

enable_language(C)
enable_language(CXX)
enable_language(Fortran)

# make the build
set(CMAKE_EXE_LINKER_FLAGS "-mcmodel=large")

set(CMAKE_EXPORT_COMPILE_COMMANDS on)

# Set up FORTRAN flags
set(CMAKE_Fortran_FLAGS "-DLINUX -Dlinux -ffixed-line-length-none -fconvert=big-endian -fno-automatic -fno-range-check -cpp -fbounds-check -g -fbacktrace -ffpe-trap=zero,overflow,underflow -std=legacy")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops")
set(CMAKE_Fortran_FLAGS_DEBUG "-ggdb -Wall -Wno-unused-variable -Wno-unused-label -Wno-unused-dummy-argument -Wno-conversion")

# set up OS-specific bits
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
  # standard Linux bits...
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_BSD_SOURCE -D_XOPEN_SOURCE=700 -D_DEFAULT_SOURCE")
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib;$ORIGIN/../opt/lib")
  set(EOS_ARCH linux)
  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
  # 64bit bits...
    set(EOS_ARCH linux64)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
    #set (CMAKE_SYSTEM_LIBRARY_PATH "/usr/lib/x86_64-linux-gnu")
    #set (CMAKE_SYSTEM_INCLUDE_PATH "/usr/include/x86_64-linux-gnu")
  endif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
# Mac bits
  set(CMAKE_MACOSX_RPATH 1)
  set(CMAKE_INSTALL_RPATH "@executable_path/../lib;@executable_path/../opt/lib")
  set(CMAKE_C_FLAGS "-DMACINTOSH ${CMAKE_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "-DMACINTOSH ${CMAKE_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "-DMACINTOSH ${CMAKE_Fortran_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "-Wl,-no_compact_unwind  -Wl,-flat_namespace ${CMAKE_EXE_LINKER_FLAGS}")
  set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-flat_namespace ${CMAKE_SHARED_LINKER_FLAGS}")
  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm64")
    set(EOS_ARCH macarm64)
  else()
    set(EOS_ARCH macintel64)
  endif()
else(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
# Miscellaneous - not sure what we are...
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib;$ORIGIN/../opt/lib")
  message("CMAKE_SYSTEM_NAME: " ${CMAKE_SYSTEM_NAME})
  message("System not recognized, using other")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

# where, oh where are your FORTRAN mod files...
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran_modules)
message("FortanModule" ${CMAKE_Fortran_MODULE_DIRECTORY})

# where we put the binaries
if(DEFINED ENV{OCSSWROOT})
  set(CMAKE_INSTALL_PREFIX $ENV{OCSSWROOT})
else()
  set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

#find_package(HDF4 REQUIRED)
#include_directories(${HDF4_INCLUDE_DIRS})
#find_package(HDFEOS REQUIRED)
#include_directories(${HDFEOS_INCLUDE_DIRS})
#set(LIBS ${CMAKE_CURRENT_SOURCE_DIR}/opt/lib/libxdr.a ${HDF4_LIBRARIES})
#set(LIBS ${HDF4_LIBRARIES})
#list(APPEND LIBS ${HDFEOS_LIBRARIES})
#message("HDF4 libraries " ${HDF4_DIR})

find_package(HDF5 REQUIRED COMPONENTS C HL Fortran Fortran_HL)
include_directories(${HDF5_INCLUDE_DIRS})
#find_package(HDFEOS5 REQUIRED)
#include_directories(${HDFEOS5_INCLUDE_DIRS})
#list(APPEND LIBS ${HDF5_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES} ${HDFEOS5_LIBRARIES})
list(APPEND LIBS ${HDF5_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})
#list(APPEND LIBS ${HDFEOS5_LIBRARIES})
message("HDF5 libraries " ${HDF5_LIBRARIES})

set (NETCDF_F77 "YES")
set (NETCDF_F90 "YES")
find_package(NetCDF REQUIRED)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_F90_LIBRARIES})
list(APPEND LIBS ${NETCDF_F77_LIBRARIES})


add_executable(oci_ua_aer
  # UV source files
  UV/src/util/ErrorHandler_class.f90
  UV/src/util/DataTypeDef.f90 
  UV/src/util/MyConstants.f90
  UV/src/util/MyMathUtil.f90
  UV/src/InterpolationModule.f90
  UV/src/util/H5Util_class.f90 
  UV/src/util/H5write_module.f90 
  UV/src/util/H5read_module.f90 
  Main/OCIUAAER_Config_Module.f90
  Main/OCIUAAER_L1BModule.f90  
  Main/Proxydata_L1BModule.f90 
  Main/Merge_land_ocean_Module.f90 
   
#  UV/src/LookupTableModule.f90  
#  UV/src/GetSurfAlbLUT_H5module.f90
#  UV/src/GetLUT_LER_AI_H5module.f90
#  UV/src/GetLUT_Miecloud_AI_H5module.f90
#  UV/src/Get_OceanLUT_H5module.f90
#  UV/src/Get_TerrainPressure_H5module.f90
#  UV/src/Get_SnowIce_module.f90 
#  UV/src/GetLUT_H5module.f90
#  UV/src/Get_ssaclimLUT2_H5module.f90
#  UV/src/Get_omacaLUT7dim_H5module.f90

  UV/src/nc4/LookupTableModule_nc4.f90  
  UV/src/nc4/GetSurfAlbLUT_H5module_nc4.f90
  UV/src/nc4/GetLUT_LER_AI_H5module_nc4.f90
  UV/src/nc4/GetLUT_Miecloud_AI_H5module_nc4.f90
  UV/src/nc4/Get_OceanLUT_H5module_nc4.f90
  UV/src/nc4/Get_TerrainPressure_H5module_nc4.f90
  UV/src/nc4/Get_SnowIce_module_nc4.f90 
  UV/src/nc4/GetLUT_H5module_nc4.f90
  UV/src/nc4/Get_ssaclimLUT2_H5module_nc4.f90
  UV/src/nc4/Get_omacaLUT7dim_H5module_nc4.f90

  DB_DT_Land/DB_Land/src/calendars.f95 
  UV/src/read_for_nuv_test.f90
  UV/src/Nearuv_alhssa_Module.f90
  UV/src/regpolymonial_predict_ssa.f90
  UV/src/NUV_AerosolModule.f90
  Main/NUV_ACAerosolModule.f90
  UV/src/NUV_package_modules.f90 
  UV/src/OCI_UV1km_DataModule.f90
  UV/src/OCI_UVAI_MieCloudModule.f90
    
  # DT Source files common to ocean and Land
  DT_Common_src/src/Interp.f90
  DT_Common_src/src/Read_gdas.f90
  DT_Common_src/src/Read_GMAO.f90 
  DT_Common_src/src/init_save_array.f90
  DT_Common_src/src/convert_rad.f90
  DT_Common_src/src/Apply_Trans_twoway.f90
  DT_Common_src/src/compute_Boxes_and_edges.f90
  DT_Common_src/src/compute_water_land.f90
  DT_Common_src/src/compute_center_GEOLOC_ANGLE.f90
  DT_Common_src/src/compute_center_Glint_Angle.f90
  DT_Common_src/src/compute_Gascorrection.f90
  DT_Common_src/src/Move_Data_Reso.f90
  DT_Common_src/src/SET_INDEX_L1b.f90
  DT_Common_src/src/open_l1b_OCI_data.f90
  DT_Common_src/src/Fillvalue_MOD04_V2.f90
  DT_Common_src/src/compute_temp_from_rad.f90       
  
  # DT Source files
  DT_Ocean/src/write_pacedt_ocean.f90
  DT_Ocean/src/write_UV.f90
  DT_Ocean/src/Handler_ocean.f90
  DT_Ocean/src/byte_set.f90
  DT_Ocean/src/compute_indx_wspeed_forLUT.f90
  DT_Ocean/src/compute_wspeed_spatial.f90
  DT_Ocean/src/Fill_Output_SDS_Ocean.f90
  DT_Ocean/src/make_ocean_cloudmask.f90
  DT_Ocean/src/Fill_Output_SDS_UV.f90
  DT_Ocean/src/Compute_Sd_for_Dust.f90
  
#  DT_Ocean/src/UV_processing.f90
#  DT_Ocean/src/Process_ocean_V2.f90
#  DT_Ocean/src/DT_package_Ocean_main.f90
#  DB_DT_Land/DT_Land/src/DT_package_Land_main.f90
#  DB_DT_Land/DT_Land/src/Process_land_V6.f90 
    
  DT_Common_src/src/nc4/dt_nc4.f90 
  DT_Common_src/src/nc4/UV_processing_nc4.f90
  DT_Common_src/src/nc4/Process_ocean_V2_nc4.f90
  DT_Common_src/src/nc4/DT_package_Ocean_main_nc4.f90
  DT_Common_src/src/nc4/DT_package_Land_main_nc4.f90
  DT_Common_src/src/nc4/Process_land_V6_nc4.f90 
  
  # DB Source files
  DB_DT_Land/DB_Land/src/pack_412.f 
  DB_DT_Land/DB_Land/src/lmpar.f 
  DB_DT_Land/DB_Land/src/fdjac2.f 
  DB_DT_Land/DB_Land/src/lmdif1.f 
  DB_DT_Land/DB_Land/src/lmdif_v2.f 
  DB_DT_Land/DB_Land/src/lmdif.f 
  DB_DT_Land/DB_Land/src/dpmpar.f 
  DB_DT_Land/DB_Land/src/enorm.f 
  DB_DT_Land/DB_Land/src/qrfac.f 
  DB_DT_Land/DB_Land/src/qrsolv.f 
  DB_DT_Land/DB_Land/src/covar.f 
  DB_DT_Land/DB_Land/src/fdjac2_v2.f 
  DB_DT_Land/DB_Land/src/atmos.f90 
  DB_DT_Land/DB_Land/src/intero.f 
  DB_DT_Land/DB_Land/src/interr.f 
  DB_DT_Land/DB_Land/src/reflec.f 
  DB_DT_Land/DB_Land/src/lodbuf.f 
  DB_DT_Land/DB_Land/src/prflec.f 
  DB_DT_Land/DB_Land/src/lodsmp.f 
  DB_DT_Land/DB_Land/src/total.f 
  DB_DT_Land/DB_Land/src/sortrx.f 
  DB_DT_Land/DB_Land/src/lmdif1_v2.f 
  DB_DT_Land/DB_Land/src/viirs_ocean_aero.f95 
  DB_DT_Land/DB_Land/src/viirs_config.f95 
  DB_DT_Land/DB_Land/src/set_qa_land.f95 
  DB_DT_Land/DB_Land/src/set_qa_ocean.f95 
  DB_DT_Land/DB_Land/src/modis.f90 
  DB_DT_Land/DB_Land/src/find_v_vegset.f 
  DB_DT_Land/DB_Land/src/find_v_viirs.f 
  DB_DT_Land/DB_Land/src/viirs_ancillary.f95 
  DB_DT_Land/DB_Land/src/screening.f95 
  
#  DB_DT_Land/DB_Land/src/DB_package_Land_main.f90
#  DB_DT_Land/DB_Land/src/read_l1b.f95 
#  DB_DT_Land/DB_Land/src/viirs_ler_luts.f95
#  DB_DT_Land/DB_Land/src/landcover.f95 
#  DB_DT_Land/DB_Land/src/get_lut_igbp_land_cover.f90 
#  DB_DT_Land/DB_Land/src/get_lut_211sfc.f90 
#  DB_DT_Land/DB_Land/src/modis_surface.f95 
#  DB_DT_Land/DB_Land/src/viirs_aerosol_luts.f95 
#  DB_DT_Land/DB_Land/src/seawifs_surface_pressure.f95 
#  DB_DT_Land/DB_Land/src/viirs_db_utils.f95 
#  DB_DT_Land/DB_Land/src/db_debug.f95 

  DB_DT_Land/DB_Land/src/nc4/DB_package_Land_main_nc4.f90
  DB_DT_Land/DB_Land/src/nc4/read_l1b_nc4.f95 
  DB_DT_Land/DB_Land/src/nc4/viirs_ler_luts_nc4.f95
  DB_DT_Land/DB_Land/src/nc4/landcover_nc4.f95
  DB_DT_Land/DB_Land/src/nc4/get_lut_igbp_land_cover_nc4.f90
  DB_DT_Land/DB_Land/src/nc4/get_lut_211sfc_nc4.f90
  DB_DT_Land/DB_Land/src/nc4/modis_surface_nc4.f95
  DB_DT_Land/DB_Land/src/nc4/viirs_aerosol_luts_nc4.f95
  DB_DT_Land/DB_Land/src/nc4/seawifs_surface_pressure_nc4.f95 
  DB_DT_Land/DB_Land/src/nc4/viirs_db_utils_nc4.f95 
#  DB_DT_Land/DB_Land/src/nc4/modis_surface_x.f95 
  
  #DT_SRCS_Land
  DB_DT_Land/DT_Land/src/compute_snow_flag.f90
  DB_DT_Land/DT_Land/src/write_pacedt_land.f90
  DB_DT_Land/DT_Land/src/compute_theshold_water_land_coastal.f90
  DB_DT_Land/DT_Land/src/make_land_cloudmask.f90
  DB_DT_Land/DT_Land/src/Handler_land.f90
  DB_DT_Land/DT_Land/src/cld_flags_forLand.f90
  DB_DT_Land/DT_Land/src/Fill_Output_SDS_land.f90
  
  #S_Write
  Write_L2_data/src/write_merged_file.f90
  Write_L2_data/src/write_1KM_merged.f90
  Write_interm_file/src/write_Output_forUVtau.f90
  Write_interm_file/src/read_Output_forUVtau.f90
  Main/OCIUAAER.f90  
)

set_target_properties(oci_ua_aer PROPERTIES
  LINK_SEARCH_START_STATIC ON
  LINK_SEARCH_END_STATIC ON
)

target_link_libraries(oci_ua_aer ${LIBS})

install (TARGETS oci_ua_aer DESTINATION bin)
