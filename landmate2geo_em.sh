#!/bin/bash

domains=$1 #("d01" "d02")
method="con"

for domain in ${domains}; do
	echo "Working on the ${domain}"

	# Extracting information for the inteprolation
	cdo selname,LANDUSEF,XLAT_M,XLONG_M geo_em.${domain}.nc LANDUSEF_WRF.nc
	ncl grid_corners.ncl

	# Interpolating LANDMATE and water_ice files to the WRf grid
	cdo remap${method},out_grid.nc LANDMATE_PFT_v1.1_Europe_0.018deg_2015.nc LANDMATE_PFT_v1.1_Europe_2015.nc
	cdo remap${method},out_grid.nc water_ice.nc ESACCI-LC-L4-LCCS_water_ice.nc

	# Compiling the fortran code for recategorization
	# NOTE: The paths to the libraries in the Makefile needs to be ajusted to your machine
	make veryclean
	make
	
	# Recategorization with the fortran program
	./WRF_LUCAS_PFTs_v3_water_ice_landmask

	# Renaming the originla file to avoud the overwriting	
	cp geo_em.${domain}.nc geo_em.${domain}_LANDMATE.nc

	# Rewriting landuse variabels 
	for varname in "LU_INDEX" "LANDUSEF" "LANDMASK"; do
		ncks -C -A -v ${varname} WRF_LUCAS_LANDMATE_PFT_v1.1_Europe_2015_${varname}_v2.nc geo_em.${domain}_LANDMATE.nc
	done

	# Delete all intermediete file
	rm WRF_LUCAS_LANDMATE_PFT_v1.1_Europe_2015_*_v2.nc
	rm ESACCI-LC-L4-LCCS_water_ice.nc
	rm LANDMATE_PFT_v1.1_Europe_2015.nc
	rm out_grid.nc LANDUSEF_WRF.nc
done


