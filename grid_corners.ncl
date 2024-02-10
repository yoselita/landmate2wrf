load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
;---Input file
    srcFileName = "LANDUSEF_WRF.nc"

;---Output file
    srcGridName = "out_grid.nc"

;----------------------------------------------------------------------
; Convert original source NCEP grid to a SCRIP convention file.
;----------------------------------------------------------------------
    src_file = addfile(srcFileName,"r")
    lat_var  = src_file->XLAT_M(0,:,:)
    lon_var  = src_file->XLONG_M(0,:,:)

    ;--- set options
    Opt                = True
    Opt@ForceOverwrite = True
    Opt@PrintTimings   = True
    Opt@Title          = "PISM Grid"
    Opt@Debug          = True

    curvilinear_to_SCRIP(srcGridName,lat_var,lon_var,Opt)

;---Clean up
    delete(Opt)

end                                             
