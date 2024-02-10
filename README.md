# List of files and description:
1. grid_corners.ncl	- ncl script that creates out_grid.nc based on your original geo_em file that is used for the interpolation with cdo
2. water_ice.nc		- nc file containing information on water and snow-ice at 300m resolution
3. WRF_LUCAS_PFTs_v3_water_ice_landmask.f90 - fortran code that performs recategorisation of the LANMATE_PFT data to the WRF grid 
4. landmate2geo_em.sh 	- bash wrapper script that performs complete job neccessary to get the final geo_em file with land cover and land mask data based on LANDMATE_PFT


# Instructions to run the script to obtain maps based on LANDMATE_PFT data
To get LANDMATE_PFT data on your WRF grid you need:
1. to place your geo_em.${domain}.nc file in your working directory, with domain="d01" (or "d02"...)
2. to copy or download LANDMATE_PFT_v1.1_Europe_0.018deg_2015.nc to the working directory
3. update Makefile with the paths to the corresponing compiler (gfortran, Intel, etc.) and netcdf libraries at you local machine
4. run ./landmate2geo_em.sh ${domian} (${domian} from point 1, e.g. d01 or d02)

After running the script, the file geo_em.${domain}_LANDMATE.nc will contain the updated LU_INDEX, LANDUSEF, and LANDMASK variables based on LANDMATE_PFT data. 

contact: 
Josipa Milovac: milovacj@unican.es 
Rita Margarida Cardoso: rmcardoso@ciencias.ulisboa.pt

