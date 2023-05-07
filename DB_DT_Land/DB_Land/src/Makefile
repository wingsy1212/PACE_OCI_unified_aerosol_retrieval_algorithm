SHELL = /bin/sh
# -fbounds-check -Wall = -Waliasing, -Wampersand, -Wconversion, -Wsurprising, -Wc-binding-type, -Wintrinsics-std, -Wno-tabs, -Wintrinsic-shadow, -Wline-truncation, -Wtarget-lifetime, -Wreal-q-constant -Wunused
#FFLAGS=-O -Waliasing -Wampersand -Wsurprising -Wintrinsics-std -Wintrinsic-shadow -Wline-truncation -Wunused -ffree-line-length-132 -ffixed-line-length-132 ${LDFLAGS}
FFLAGS=-O -ffixed-line-length-132 ${LDFLAGS}
LFLAGS = -lmfhdf -ldf -lhdf5hl_fortran -lhdf5_hl -lnetcdff -lnetcdf -lhdf5_fortran -lhdf5 -ljpeg -lz -lgrib_api -lgrib_api_f90

.PHONY: clean

deep_blue : sort.o median.o pack_412.o lmdif1.o lmdif.o dpmpar.o enorm.o fdjac2.o lmpar.o qrfac.o qrsolv.o covar.o fdjac2_v2.o lmdif1_v2.o lmdif_v2.o atmos.o calendars.o intero.o interr.o reflec.o viirs_ler_luts.o lodbuf.o prflec.o lodsmp.o total.o sortrx.o landcover.o viirs_obpg_corrections.o viirs_aerosol_luts.o viirs_ocean_aero.o seawifs_surface_pressure.o modis_surface.o viirs_config.o viirs_db_utils.o read_data.o set_qa_land.o set_qa_ocean.o read_l1b.o modis.o get_lut_211sfc.o get_lut_igbp_land_cover.o db_debug.o find_v_vegset.o find_v_viirs.o viirs_ancillary.o screening.o  deep_blue.f90
	${FC} ${FFLAGS} ${LDFLAGS} $^ -o $@ ${LFLAGS}
read_l1b.o : read_l1b.f95
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}			
test_nc4 : test_nc4.f95
	${FC} ${FFLAGS} ${LDFLAGS} $^ -o $@ ${LFLAGS}
modis_surface.o : sortrx.o modis_surface.f95
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
viirs_ancillary.o : viirs_ancillary.f95
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
db_debug.o : db_debug.f95
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
viirs_obpg_corrections.o : viirs_obpg_corrections.f95
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
viirs_db_utils.o : viirs_db_utils.f95
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
viirs_config.o : viirs_config.f95
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
landcover.o : landcover.f95
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
seawifs_surface_pressure.o : seawifs_surface_pressure.f95
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
modis.o : modis.f90
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
read_data.o : read_data.f
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
set_qa_land.o : set_qa_land.f95
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
set_qa_ocean.o : set_qa_ocean.f95
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
pack_412.o : pack_412.f
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
viirs_aerosol_luts.o : viirs_aerosol_luts.f95
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
viirs_ocean_aero.o : median.o lmdif1.o viirs_ocean_aero.f95
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
viirs_ler_luts.o : viirs_ler_luts.f95
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
find_v_viirs.o : calendars.o find_v_viirs.f
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
find_v_vegset.o : find_v_vegset.f
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
get_lut_igbp_land_cover.o : get_lut_igbp_land_cover.f90
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
get_lut_211sfc.o : get_lut_211sfc.f90
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
sortrx.o : sortrx.f
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
calendars.o : calendars.f95
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
total.o : total.f
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
lodsmp.o : lodsmp.f
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
prflec.o : prflec.f
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
lodbuf.o : lodbuf.f
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
reflec.o : interr.o intero.o reflec.f
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
interr.o : interr.f
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
intero.o : intero.f
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
atmos.o : atmos.f90
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
screening.o : screening.f95
	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}		
lmdif1.o : dpmpar.o enorm.o fdjac2.o lmpar.o qrfac.o qrsolv.o lmdif.o lmdif1.f
	${FC} ${FFLAGS} -c $^
dpmpar.o : dpmpar.f
	${FC} ${FFLAGS} -c $^
enorm.o : enorm.f
	${FC} ${FFLAGS} -c $^
fdjac2.o : fdjac2.f
	${FC} ${FFLAGS} -c $^
lmpar.o : lmpar.f
	${FC} ${FFLAGS} -c $^
qrfac.o : qrfac.f
	${FC} ${FFLAGS} -c $^
qrsolv.o : qrsolv.f
	${FC} ${FFLAGS} -c $^
lmdif.o : lmdif.f
	${FC} ${FFLAGS} -c $^
lmdif1_v2.o : dpmpar.o enorm.o fdjac2_v2.o lmpar.o qrfac.o qrsolv.o lmdif_v2.o lmdif1_v2.f
	${FC} ${FFLAGS} -c $^
fdjac2_v2.o : fdjac2_v2.f
	${FC} ${FFLAGS} -c $^   
lmdif_v2.o : lmdif_v2.f
	${FC} ${FFLAGS} -c $^
covar.o : covar.f
	${FC} ${FFLAGS} -c $^      	
median.o : sort.o median.f
	${FC} ${LDFLAGS} -c $^
sort.o : sort.f
	${FC} ${LDFLAGS} -c $^
create_aerosol_model_input : calendars.o create_aerosol_model_input.f95
	${FC} ${FFLAGS} ${LDFLAGS} $^ -o $@ ${LFLAGS}
create_brdfbase_input : calendars.o create_brdfbase_input.f95
	${FC} ${FFLAGS} ${LDFLAGS} $^ -o $@ ${LFLAGS}
create_geozone_input : calendars.o create_geozone_input.f95
	${FC} ${FFLAGS} ${LDFLAGS} $^ -o $@ ${LFLAGS}


l2_merge : calendars.o viirs_config.o viirs_db_utils.o l2_merge.f90
	${FC} ${FFLAGS} ${LDFLAGS} $^ -o $@ ${LFLAGS}
# viirs_db_utils.o : viirs_db_utils.f95
# 	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
# viirs_config.o : viirs_config.f95
# 	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}
# calendars.o : calendars.f95
# 	${FC} -c ${FFLAGS} ${LDFLAGS} $^ ${LFLAGS}	
	
clean :
	rm sort.o median.o pack_412.o lmdif1.o lmdif.o dpmpar.o enorm.o fdjac2.o lmpar.o qrfac.o \
	qrsolv.o atmos.o calendars.o intero.o interr.o reflec.o viirs_ler_luts.o lodbuf.o prflec.o \
	lodsmp.o total.o sortrx.o landcover.o viirs_obpg_corrections.o viirs_aerosol_luts.o \
	viirs_ocean_aero.o seawifs_surface_pressure.o modis_surface.o viirs_config.o viirs_db_utils.o \
	read_data.o set_qa_land.o set_qa_ocean.o read_l1b.o modis.o get_lut_211sfc.o get_lut_igbp_land_cover.o \
	db_debug.o find_v_vegset.o find_v_viirs.o viirs_ancillary.o screening.o deep_blue
