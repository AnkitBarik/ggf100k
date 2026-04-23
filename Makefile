FC = gfortran
FFLAGS = -O3 -fPIC -fcheck=all -Wall -Wextra -Werror
HDF5_LIBS = -L/usr/lib/x86_64-linux-gnu/ -lhdf5_serial_fortran
HDF5_INCLUDES = -I/usr/include/hdf5/serial/

SRCS = getCoeffs.f90

getCoeffs.exe : $(SRCS)
	$(FC) $(FFLAGS) $(SRCS) -o $@ $(HDF5_LIBS) $(HDF5_INCLUDES)

all: getCoeffs.exe
clean:
	rm -f getCoeffs.exe *.mod

.PHONY: all clean