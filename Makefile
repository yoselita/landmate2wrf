#Ubuntu Desktop, -g or -O2
FC = /oceano/gmeteo/users/milovacj/miniconda3/envs/netcdf/bin/gfortran
FCFLAGS = -O2
#FCFLAGS += -Wall
FCFLAGS += -ffree-line-length-none
FCFLAGS += -Wno-tabs
FCFLAGS += -I/oceano/gmeteo/users/milovacj/miniconda3/envs/netcdf/include/
LDFLAGS = -L/oceano/gmeteo/users/milovacj/miniconda3/envs/netcdf/lib -lnetcdff -lnetcdf


PROGRAMS = WRF_LUCAS_PFTs_v3_water_ice_landmask 

all: $(PROGRAMS)

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD *_genmod.f90 tmpfile* log

veryclean: clean
	rm -rf *~ $(PROGRAMS)

